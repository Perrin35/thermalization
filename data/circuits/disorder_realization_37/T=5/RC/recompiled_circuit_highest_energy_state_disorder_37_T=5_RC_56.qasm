OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3647519) q[0];
sx q[0];
rz(-2.2651894) q[0];
sx q[0];
rz(-2.8887698) q[0];
rz(-1.9934935) q[1];
sx q[1];
rz(-0.91528457) q[1];
sx q[1];
rz(0.80438703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1348477) q[0];
sx q[0];
rz(-0.50323707) q[0];
sx q[0];
rz(-0.87849599) q[0];
rz(-pi) q[1];
rz(-2.7115066) q[2];
sx q[2];
rz(-1.7835435) q[2];
sx q[2];
rz(1.9274005) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6614395) q[1];
sx q[1];
rz(-0.51750187) q[1];
sx q[1];
rz(-0.75019849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31320235) q[3];
sx q[3];
rz(-2.7814908) q[3];
sx q[3];
rz(-1.2822156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.048451) q[2];
sx q[2];
rz(-2.1790049) q[2];
sx q[2];
rz(0.87240458) q[2];
rz(0.33556542) q[3];
sx q[3];
rz(-1.395697) q[3];
sx q[3];
rz(-0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5035079) q[0];
sx q[0];
rz(-2.8758949) q[0];
sx q[0];
rz(0.22915325) q[0];
rz(2.9511662) q[1];
sx q[1];
rz(-1.8016022) q[1];
sx q[1];
rz(0.71358877) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021727) q[0];
sx q[0];
rz(-1.3029126) q[0];
sx q[0];
rz(-2.2557206) q[0];
rz(-pi) q[1];
rz(2.607478) q[2];
sx q[2];
rz(-0.91098173) q[2];
sx q[2];
rz(0.5420064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4137607) q[1];
sx q[1];
rz(-1.276386) q[1];
sx q[1];
rz(2.4319699) q[1];
x q[2];
rz(0.76477625) q[3];
sx q[3];
rz(-1.4524609) q[3];
sx q[3];
rz(1.8272994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4552292) q[2];
sx q[2];
rz(-0.31193048) q[2];
sx q[2];
rz(2.6658106) q[2];
rz(2.0717715) q[3];
sx q[3];
rz(-2.1333623) q[3];
sx q[3];
rz(-0.76655918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0723202) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(1.2323761) q[0];
rz(2.6909289) q[1];
sx q[1];
rz(-1.181299) q[1];
sx q[1];
rz(-0.88952363) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7018676) q[0];
sx q[0];
rz(-2.151818) q[0];
sx q[0];
rz(2.0608611) q[0];
x q[1];
rz(0.81368229) q[2];
sx q[2];
rz(-0.70975862) q[2];
sx q[2];
rz(-0.26772945) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0625713) q[1];
sx q[1];
rz(-1.3600742) q[1];
sx q[1];
rz(0.91075588) q[1];
rz(1.9185478) q[3];
sx q[3];
rz(-1.2158211) q[3];
sx q[3];
rz(2.6048911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.687261) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.8827776) q[2];
rz(-1.2339633) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(-0.0080571938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71037978) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(1.4696962) q[0];
rz(1.1471033) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(-2.240644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44997901) q[0];
sx q[0];
rz(-1.8834582) q[0];
sx q[0];
rz(-1.6603966) q[0];
x q[1];
rz(0.9587165) q[2];
sx q[2];
rz(-1.0462073) q[2];
sx q[2];
rz(0.9768578) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4576396) q[1];
sx q[1];
rz(-1.4767273) q[1];
sx q[1];
rz(2.4165618) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4114831) q[3];
sx q[3];
rz(-0.84419227) q[3];
sx q[3];
rz(-0.58578736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49742302) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(-0.72745848) q[2];
rz(-0.71508956) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(2.6434744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6775976) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(2.4053307) q[0];
rz(0.62034208) q[1];
sx q[1];
rz(-0.54519975) q[1];
sx q[1];
rz(-1.6166519) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10337457) q[0];
sx q[0];
rz(-3.0854825) q[0];
sx q[0];
rz(-0.17608626) q[0];
rz(-pi) q[1];
rz(-1.9061766) q[2];
sx q[2];
rz(-1.35193) q[2];
sx q[2];
rz(0.21974213) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1279631) q[1];
sx q[1];
rz(-1.9444325) q[1];
sx q[1];
rz(-2.29261) q[1];
x q[2];
rz(-0.81186632) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(-0.52385274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6058495) q[2];
sx q[2];
rz(-0.77631408) q[2];
sx q[2];
rz(2.1898451) q[2];
rz(-2.7239299) q[3];
sx q[3];
rz(-0.2499191) q[3];
sx q[3];
rz(-1.1550268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8010537) q[0];
sx q[0];
rz(-1.2507573) q[0];
sx q[0];
rz(-2.387555) q[0];
rz(0.89318371) q[1];
sx q[1];
rz(-2.1008396) q[1];
sx q[1];
rz(-2.975614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245859) q[0];
sx q[0];
rz(-2.4501194) q[0];
sx q[0];
rz(2.9777479) q[0];
rz(1.7506261) q[2];
sx q[2];
rz(-1.0713312) q[2];
sx q[2];
rz(1.3701554) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2365723) q[1];
sx q[1];
rz(-0.81099866) q[1];
sx q[1];
rz(-0.6375957) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2955722) q[3];
sx q[3];
rz(-1.4365734) q[3];
sx q[3];
rz(-1.272066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37504998) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(-0.005391187) q[2];
rz(-3.052875) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(0.79571342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8832815) q[0];
sx q[0];
rz(-1.9323876) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(2.2329277) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(-0.044513449) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9712898) q[0];
sx q[0];
rz(-1.4409541) q[0];
sx q[0];
rz(1.4128039) q[0];
x q[1];
rz(-0.037864121) q[2];
sx q[2];
rz(-1.711039) q[2];
sx q[2];
rz(-0.87272206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9981838) q[1];
sx q[1];
rz(-2.1381731) q[1];
sx q[1];
rz(-2.8860249) q[1];
rz(-0.022734108) q[3];
sx q[3];
rz(-1.4511381) q[3];
sx q[3];
rz(-2.5496063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8160416) q[2];
sx q[2];
rz(-0.52246919) q[2];
sx q[2];
rz(2.6204056) q[2];
rz(2.9442287) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(-0.50281966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40081438) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(-0.25522301) q[0];
rz(2.9995645) q[1];
sx q[1];
rz(-1.4165001) q[1];
sx q[1];
rz(-2.7878888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86310327) q[0];
sx q[0];
rz(-1.2444087) q[0];
sx q[0];
rz(1.9410237) q[0];
rz(-pi) q[1];
rz(-0.94714099) q[2];
sx q[2];
rz(-1.9883687) q[2];
sx q[2];
rz(3.0018978) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6908011) q[1];
sx q[1];
rz(-2.004068) q[1];
sx q[1];
rz(2.7669897) q[1];
rz(1.0068302) q[3];
sx q[3];
rz(-2.2368397) q[3];
sx q[3];
rz(0.49426916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8346617) q[2];
sx q[2];
rz(-0.31535172) q[2];
sx q[2];
rz(1.5241415) q[2];
rz(-2.2751685) q[3];
sx q[3];
rz(-1.4640936) q[3];
sx q[3];
rz(-2.5518937) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50188142) q[0];
sx q[0];
rz(-0.15637936) q[0];
sx q[0];
rz(1.7891275) q[0];
rz(-1.2347219) q[1];
sx q[1];
rz(-0.23209485) q[1];
sx q[1];
rz(2.629705) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3054661) q[0];
sx q[0];
rz(-0.16634596) q[0];
sx q[0];
rz(1.8977099) q[0];
x q[1];
rz(-2.9779766) q[2];
sx q[2];
rz(-2.3694042) q[2];
sx q[2];
rz(2.8761187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3583957) q[1];
sx q[1];
rz(-1.8640717) q[1];
sx q[1];
rz(2.6986928) q[1];
x q[2];
rz(-2.0824349) q[3];
sx q[3];
rz(-2.4928204) q[3];
sx q[3];
rz(-1.2973451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1341683) q[2];
sx q[2];
rz(-1.8210501) q[2];
sx q[2];
rz(-0.39719886) q[2];
rz(-1.0154137) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(0.5526244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7427202) q[0];
sx q[0];
rz(-1.2074559) q[0];
sx q[0];
rz(2.7040828) q[0];
rz(2.4028026) q[1];
sx q[1];
rz(-2.3414325) q[1];
sx q[1];
rz(-1.4791666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.096232) q[0];
sx q[0];
rz(-2.6511565) q[0];
sx q[0];
rz(-0.11842056) q[0];
rz(-3.0257341) q[2];
sx q[2];
rz(-1.2136974) q[2];
sx q[2];
rz(3.0810205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38592713) q[1];
sx q[1];
rz(-3.0184221) q[1];
sx q[1];
rz(2.6120899) q[1];
rz(-pi) q[2];
rz(-0.78082943) q[3];
sx q[3];
rz(-2.1159647) q[3];
sx q[3];
rz(-0.050267537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9208357) q[2];
sx q[2];
rz(-1.7831384) q[2];
sx q[2];
rz(2.9662568) q[2];
rz(-1.0026503) q[3];
sx q[3];
rz(-0.23313871) q[3];
sx q[3];
rz(1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6709082) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(-2.9910174) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(0.038240959) q[2];
sx q[2];
rz(-0.19842166) q[2];
sx q[2];
rz(-1.7501065) q[2];
rz(-0.74674947) q[3];
sx q[3];
rz(-2.8058553) q[3];
sx q[3];
rz(-0.6466858) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
