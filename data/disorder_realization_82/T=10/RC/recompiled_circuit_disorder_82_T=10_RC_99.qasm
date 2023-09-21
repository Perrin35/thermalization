OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5582433) q[0];
sx q[0];
rz(-1.3238751) q[0];
sx q[0];
rz(-0.14351828) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7202397) q[2];
sx q[2];
rz(-1.7755277) q[2];
sx q[2];
rz(2.8533964) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0496088) q[1];
sx q[1];
rz(-0.99124747) q[1];
sx q[1];
rz(-2.0642573) q[1];
x q[2];
rz(0.7198556) q[3];
sx q[3];
rz(-2.8350283) q[3];
sx q[3];
rz(2.7048064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(-2.4025829) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-0.48746902) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2330403) q[0];
sx q[0];
rz(-1.4899583) q[0];
sx q[0];
rz(-3.0768865) q[0];
x q[1];
rz(2.3357453) q[2];
sx q[2];
rz(-0.57493756) q[2];
sx q[2];
rz(0.88428674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4844955) q[1];
sx q[1];
rz(-2.8946981) q[1];
sx q[1];
rz(-2.2730278) q[1];
rz(-pi) q[2];
rz(-3.0549166) q[3];
sx q[3];
rz(-1.254734) q[3];
sx q[3];
rz(2.812127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(-3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19668002) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(2.9512067) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-2.9188459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6762786) q[0];
sx q[0];
rz(-1.197581) q[0];
sx q[0];
rz(-2.7200384) q[0];
rz(-2.9163755) q[2];
sx q[2];
rz(-1.2080964) q[2];
sx q[2];
rz(-0.97937102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0157156) q[1];
sx q[1];
rz(-1.5029969) q[1];
sx q[1];
rz(-2.5281639) q[1];
rz(0.95007105) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(1.0328968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(-0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.7383204) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5266787) q[0];
sx q[0];
rz(-2.3290312) q[0];
sx q[0];
rz(-2.148669) q[0];
x q[1];
rz(1.4843066) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(0.53777018) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4184119) q[1];
sx q[1];
rz(-0.81058093) q[1];
sx q[1];
rz(-2.4852738) q[1];
x q[2];
rz(-1.7761049) q[3];
sx q[3];
rz(-1.218154) q[3];
sx q[3];
rz(-2.7384788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.197864) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.8614004) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.7368332) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(2.6884902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26950726) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(2.0399658) q[0];
x q[1];
rz(-0.034899072) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(-2.3226483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0616152) q[1];
sx q[1];
rz(-1.6528168) q[1];
sx q[1];
rz(-1.3498989) q[1];
rz(2.1987678) q[3];
sx q[3];
rz(-0.70126611) q[3];
sx q[3];
rz(3.0758465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0129464) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033427514) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(-0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(1.5135117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0673163) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(-0.61074722) q[0];
rz(-1.6983301) q[2];
sx q[2];
rz(-2.6886534) q[2];
sx q[2];
rz(-0.73308257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7442419) q[1];
sx q[1];
rz(-2.2480818) q[1];
sx q[1];
rz(2.1703297) q[1];
rz(-1.120818) q[3];
sx q[3];
rz(-1.9213772) q[3];
sx q[3];
rz(1.9812802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(-0.13051662) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(1.827084) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.1869173) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884739) q[0];
sx q[0];
rz(-1.0736199) q[0];
sx q[0];
rz(2.8209646) q[0];
rz(-pi) q[1];
x q[1];
rz(1.600012) q[2];
sx q[2];
rz(-1.511682) q[2];
sx q[2];
rz(2.7616449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4759051) q[1];
sx q[1];
rz(-2.396282) q[1];
sx q[1];
rz(1.7687294) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8625453) q[3];
sx q[3];
rz(-2.1803133) q[3];
sx q[3];
rz(-2.9677344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(-1.5073744) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54661575) q[0];
sx q[0];
rz(-0.82406509) q[0];
sx q[0];
rz(2.0440408) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0975935) q[2];
sx q[2];
rz(-1.6128522) q[2];
sx q[2];
rz(0.85925697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.028138782) q[1];
sx q[1];
rz(-2.8061562) q[1];
sx q[1];
rz(2.0466652) q[1];
x q[2];
rz(-0.52328531) q[3];
sx q[3];
rz(-2.0998294) q[3];
sx q[3];
rz(-0.39213359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(-2.3146546) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6475911) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(2.0064328) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.1816918) q[0];
sx q[0];
rz(-1.7745716) q[0];
x q[1];
rz(-1.4831545) q[2];
sx q[2];
rz(-0.32155514) q[2];
sx q[2];
rz(1.1631539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70442048) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(1.0874332) q[1];
x q[2];
rz(-1.5376066) q[3];
sx q[3];
rz(-0.6233218) q[3];
sx q[3];
rz(0.25191307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(-0.19616729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575532) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(-2.2474225) q[0];
rz(-pi) q[1];
rz(0.050373366) q[2];
sx q[2];
rz(-1.435624) q[2];
sx q[2];
rz(-0.86063517) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4539459) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(-0.83666283) q[1];
x q[2];
rz(0.85828652) q[3];
sx q[3];
rz(-0.86601102) q[3];
sx q[3];
rz(1.5376877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(-1.1817415) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186196) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(1.3148057) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(-2.9076004) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(1.8004988) q[3];
sx q[3];
rz(-1.5256186) q[3];
sx q[3];
rz(-1.516173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
