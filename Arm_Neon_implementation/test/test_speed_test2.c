#include <stdint.h>
#include "../sign.h"
#include "../sign_ntt_c.h"
#include "../sign_ntt_neon.h"
#include "../sign_spm_neon.h"
#include "../sign_hybrid_neon.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../params.h"
#include "cpucycles.h"
#include "speed_print.h"

#include <linux/perf_event.h>
#include <linux/hw_breakpoint.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>

#include <time.h>

uint64_t getCpuCycles(void){
  struct timespec time;

  clock_gettime(CLOCK_REALTIME, &time);
  return (uint64_t)(time.tv_sec*1e9 + time.tv_nsec);
}

struct read_format {
    uint64_t nr;
    struct {
        uint64_t value;
        uint64_t id;
    } values[];
};

int create_hardware_perf(int grp_fd, enum perf_hw_id hw_ids, uint64_t *ioc_id)
{
    if(PERF_COUNT_HW_MAX <= hw_ids || hw_ids < 0) {
        printf("Unsupport enum perf_hw_id.\n");
        return -1;
    }
    
    struct perf_event_attr pea;
    
    memset(&pea, 0, sizeof(struct perf_event_attr));
    pea.type = PERF_TYPE_HARDWARE;
    pea.size = sizeof(struct perf_event_attr);
    pea.config = hw_ids;
    pea.disabled = 1;
    pea.exclude_kernel = 1;
    pea.exclude_hv = 1;
    //pea.read_format = PERF_FORMAT_ID;
    int fd = syscall(__NR_perf_event_open, &pea, 0, -1, grp_fd>2?grp_fd:-1, 0);
    //ioctl(fd, PERF_EVENT_IOC_ID, ioc_id);

    return fd;
}

#define NTESTS 1000

uint64_t t[NTESTS];

int main(void)
{
  unsigned int i;
  size_t smlen;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sm[CRYPTO_BYTES + CRHBYTES];

  int group_fd;
  uint64_t id1;
  //char buf[4096];
  //struct read_format* rf = (struct read_format*) buf;
  group_fd = create_hardware_perf(-1, PERF_COUNT_HW_CPU_CYCLES, &id1);

  ioctl(group_fd, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    crypto_sign_keypair(pk, sk);
  }
  print_results("Keypair:", t, NTESTS);

  // 1.ntt c implementation
  printf("---ntt_c---\n");
  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    ntt_c_crypto_sign(sm, &smlen, sm, CRHBYTES, sk);
  }
  print_results("Sign:", t, NTESTS);

  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    ntt_c_crypto_sign_verify(sm, CRYPTO_BYTES, sm, CRHBYTES, pk);
  }
  print_results("Verify:", t, NTESTS);

  // 2.ntt neon
  printf("---ntt_neon---\n");
  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    ntt_neon_crypto_sign(sm, &smlen, sm, CRHBYTES, sk);
  }
  print_results("Sign:", t, NTESTS);

  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    ntt_neon_crypto_sign_verify(sm, CRYPTO_BYTES, sm, CRHBYTES, pk);
  }
  print_results("Verify:", t, NTESTS);

  // 3.parallel small polynomial multiplication neon, spm_neon
  printf("---spm_neon---\n");
  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    spm_neon_crypto_sign(sm, &smlen, sm, CRHBYTES, sk);
  }
  print_results("Sign:", t, NTESTS);

  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    spm_neon_crypto_sign_verify(sm, CRYPTO_BYTES, sm, CRHBYTES, pk);
  }
  print_results("Verify:", t, NTESTS);

  // 4.hybrid neon, hybrid_neon
  printf("---hybrid_neon: spm_neon & ntt_neon---\n");
  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    hybrid_neon_crypto_sign(sm, &smlen, sm, CRHBYTES, sk);
  }
  print_results("Sign:", t, NTESTS);

  ioctl(group_fd, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  for(i = 0; i < NTESTS; ++i) {
    read(group_fd, &t[i], sizeof(t[i]));
    hybrid_neon_crypto_sign_verify(sm, CRYPTO_BYTES, sm, CRHBYTES, pk);
  }
  print_results("Verify:", t, NTESTS);

  ioctl(group_fd, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
  close(group_fd);
  return 0;
}
