OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(-1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(-3.1187305) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6883144) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(-0.90318371) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1360264) q[2];
sx q[2];
rz(-2.044894) q[2];
sx q[2];
rz(-2.9413162) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3373128) q[1];
sx q[1];
rz(-0.37706456) q[1];
sx q[1];
rz(-1.1537329) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2223577) q[3];
sx q[3];
rz(-2.493353) q[3];
sx q[3];
rz(-0.1227619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(0.76618761) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(3.0564953) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3413977) q[0];
sx q[0];
rz(-1.5543823) q[0];
sx q[0];
rz(-0.62231681) q[0];
rz(-pi) q[1];
rz(-0.9688103) q[2];
sx q[2];
rz(-1.8245056) q[2];
sx q[2];
rz(-1.837455) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47463402) q[1];
sx q[1];
rz(-0.14161319) q[1];
sx q[1];
rz(0.9160362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3195023) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(3.086123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.3844301) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2565838) q[0];
sx q[0];
rz(-0.48261595) q[0];
sx q[0];
rz(-2.9628739) q[0];
rz(-pi) q[1];
rz(2.7483814) q[2];
sx q[2];
rz(-2.0438497) q[2];
sx q[2];
rz(-0.87393239) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78450173) q[1];
sx q[1];
rz(-2.38584) q[1];
sx q[1];
rz(-3.109039) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.609385) q[3];
sx q[3];
rz(-2.0720707) q[3];
sx q[3];
rz(-1.7621079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79065381) q[0];
sx q[0];
rz(-1.5293984) q[0];
sx q[0];
rz(-0.008965094) q[0];
x q[1];
rz(-1.3992594) q[2];
sx q[2];
rz(-2.9651387) q[2];
sx q[2];
rz(-1.2982969) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85992766) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(1.330919) q[1];
rz(-pi) q[2];
x q[2];
rz(1.585235) q[3];
sx q[3];
rz(-2.0824021) q[3];
sx q[3];
rz(-2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8308782) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(0.80835289) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-2.2461058) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(-0.53363824) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424608) q[0];
sx q[0];
rz(-1.2958382) q[0];
sx q[0];
rz(-1.2752227) q[0];
rz(-pi) q[1];
rz(-1.7970423) q[2];
sx q[2];
rz(-2.3257253) q[2];
sx q[2];
rz(-2.0362542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4662019) q[1];
sx q[1];
rz(-1.7168105) q[1];
sx q[1];
rz(-0.22202613) q[1];
rz(-2.93612) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(1.458414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50826532) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-2.1369381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6509103) q[0];
sx q[0];
rz(-2.6255529) q[0];
sx q[0];
rz(0.74923058) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(-2.4900988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55014729) q[1];
sx q[1];
rz(-2.1131556) q[1];
sx q[1];
rz(1.28246) q[1];
rz(-pi) q[2];
rz(-1.5980914) q[3];
sx q[3];
rz(-0.41397646) q[3];
sx q[3];
rz(-0.15347813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7375609) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(2.693434) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(2.4694494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7710966) q[0];
sx q[0];
rz(-1.578195) q[0];
sx q[0];
rz(1.4047983) q[0];
x q[1];
rz(2.7156098) q[2];
sx q[2];
rz(-1.9773215) q[2];
sx q[2];
rz(-2.4533518) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8183648) q[1];
sx q[1];
rz(-2.0548901) q[1];
sx q[1];
rz(1.1546043) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4268731) q[3];
sx q[3];
rz(-0.87464273) q[3];
sx q[3];
rz(-1.5211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(0.80668443) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8950302) q[0];
sx q[0];
rz(-0.44996214) q[0];
sx q[0];
rz(1.8438086) q[0];
x q[1];
rz(-2.9358747) q[2];
sx q[2];
rz(-1.4652025) q[2];
sx q[2];
rz(0.92668698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.403277) q[1];
sx q[1];
rz(-0.74790955) q[1];
sx q[1];
rz(-0.43055375) q[1];
rz(0.032012352) q[3];
sx q[3];
rz(-1.8714928) q[3];
sx q[3];
rz(-0.34919448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(0.99964833) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(-2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-2.8009169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1239615) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(2.7244363) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0057783) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(0.68067683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68742467) q[1];
sx q[1];
rz(-1.3892281) q[1];
sx q[1];
rz(1.0086235) q[1];
x q[2];
rz(-2.4386028) q[3];
sx q[3];
rz(-2.8841416) q[3];
sx q[3];
rz(1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-0.27858946) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(0.07671193) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0893223) q[0];
sx q[0];
rz(-2.8025586) q[0];
sx q[0];
rz(1.9498755) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99761988) q[2];
sx q[2];
rz(-1.7152889) q[2];
sx q[2];
rz(1.0695374) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50347603) q[1];
sx q[1];
rz(-1.5917115) q[1];
sx q[1];
rz(2.7760844) q[1];
rz(-pi) q[2];
rz(2.6807908) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(-2.868696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(-0.43141836) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(0.063522696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(2.5554399) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(-1.099844) q[2];
sx q[2];
rz(-1.8269314) q[2];
sx q[2];
rz(-1.6614428) q[2];
rz(2.7145731) q[3];
sx q[3];
rz(-1.2357161) q[3];
sx q[3];
rz(-1.2824035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
