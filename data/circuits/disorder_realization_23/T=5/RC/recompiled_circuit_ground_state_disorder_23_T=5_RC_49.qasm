OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5705559) q[0];
sx q[0];
rz(-0.67222995) q[0];
sx q[0];
rz(-1.6663405) q[0];
rz(-2.1885459) q[1];
sx q[1];
rz(-0.16521984) q[1];
sx q[1];
rz(1.8728949) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8333913) q[0];
sx q[0];
rz(-1.1723639) q[0];
sx q[0];
rz(1.3223668) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6495291) q[2];
sx q[2];
rz(-2.1644219) q[2];
sx q[2];
rz(0.77897859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70763146) q[1];
sx q[1];
rz(-0.58506706) q[1];
sx q[1];
rz(-1.8088887) q[1];
rz(-pi) q[2];
rz(2.3846486) q[3];
sx q[3];
rz(-1.401541) q[3];
sx q[3];
rz(1.7361803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98422009) q[2];
sx q[2];
rz(-2.7998388) q[2];
sx q[2];
rz(-2.3483707) q[2];
rz(0.67304099) q[3];
sx q[3];
rz(-1.0854191) q[3];
sx q[3];
rz(-1.2696772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9840045) q[0];
sx q[0];
rz(-2.3389811) q[0];
sx q[0];
rz(-2.5322835) q[0];
rz(0.20978236) q[1];
sx q[1];
rz(-2.3502626) q[1];
sx q[1];
rz(-0.9598859) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6580556) q[0];
sx q[0];
rz(-1.6153045) q[0];
sx q[0];
rz(-3.0625212) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069931937) q[2];
sx q[2];
rz(-2.0581823) q[2];
sx q[2];
rz(-2.5068381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34997019) q[1];
sx q[1];
rz(-0.73957878) q[1];
sx q[1];
rz(1.6639568) q[1];
x q[2];
rz(-2.2678492) q[3];
sx q[3];
rz(-2.5140719) q[3];
sx q[3];
rz(-2.2604159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7858872) q[2];
sx q[2];
rz(-0.83928078) q[2];
sx q[2];
rz(0.4915702) q[2];
rz(3.135904) q[3];
sx q[3];
rz(-1.0892884) q[3];
sx q[3];
rz(-0.51386851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095057644) q[0];
sx q[0];
rz(-2.6061366) q[0];
sx q[0];
rz(2.3989578) q[0];
rz(0.62478089) q[1];
sx q[1];
rz(-0.94015986) q[1];
sx q[1];
rz(2.0754441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91573533) q[0];
sx q[0];
rz(-1.7122147) q[0];
sx q[0];
rz(-1.3456324) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0451308) q[2];
sx q[2];
rz(-2.774775) q[2];
sx q[2];
rz(-2.9638673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9650363) q[1];
sx q[1];
rz(-1.6966469) q[1];
sx q[1];
rz(-1.9767889) q[1];
rz(-0.2188628) q[3];
sx q[3];
rz(-1.0893679) q[3];
sx q[3];
rz(2.4784513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8891958) q[2];
sx q[2];
rz(-2.9434581) q[2];
sx q[2];
rz(-0.55257094) q[2];
rz(2.6342454) q[3];
sx q[3];
rz(-2.0523968) q[3];
sx q[3];
rz(-0.56940091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4813389) q[0];
sx q[0];
rz(-0.57732552) q[0];
sx q[0];
rz(1.8950155) q[0];
rz(1.735911) q[1];
sx q[1];
rz(-0.77189267) q[1];
sx q[1];
rz(-0.71800047) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3221949) q[0];
sx q[0];
rz(-0.85032636) q[0];
sx q[0];
rz(2.2706991) q[0];
x q[1];
rz(2.1338855) q[2];
sx q[2];
rz(-1.3674842) q[2];
sx q[2];
rz(-0.85693923) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4990684) q[1];
sx q[1];
rz(-0.82083265) q[1];
sx q[1];
rz(-1.3715416) q[1];
x q[2];
rz(-0.67881949) q[3];
sx q[3];
rz(-2.9390237) q[3];
sx q[3];
rz(0.55965078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2672853) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(-0.70449746) q[2];
rz(0.38257515) q[3];
sx q[3];
rz(-1.4380598) q[3];
sx q[3];
rz(-2.2883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48069561) q[0];
sx q[0];
rz(-0.040204164) q[0];
sx q[0];
rz(0.30886343) q[0];
rz(0.94912306) q[1];
sx q[1];
rz(-0.32052761) q[1];
sx q[1];
rz(2.5905051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8723485) q[0];
sx q[0];
rz(-1.1820894) q[0];
sx q[0];
rz(1.7354911) q[0];
x q[1];
rz(3.0727714) q[2];
sx q[2];
rz(-2.2274979) q[2];
sx q[2];
rz(-2.7842885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84173735) q[1];
sx q[1];
rz(-1.7256712) q[1];
sx q[1];
rz(-2.9938404) q[1];
rz(-pi) q[2];
rz(2.9948009) q[3];
sx q[3];
rz(-2.3596968) q[3];
sx q[3];
rz(1.1679389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.71517313) q[2];
sx q[2];
rz(-0.74843633) q[2];
sx q[2];
rz(-0.61881649) q[2];
rz(2.5185781) q[3];
sx q[3];
rz(-2.2996733) q[3];
sx q[3];
rz(-0.4272517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039336786) q[0];
sx q[0];
rz(-0.66348851) q[0];
sx q[0];
rz(-0.46184194) q[0];
rz(0.49679187) q[1];
sx q[1];
rz(-1.8180314) q[1];
sx q[1];
rz(-1.7740446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753938) q[0];
sx q[0];
rz(-2.4552422) q[0];
sx q[0];
rz(-1.4453056) q[0];
rz(-pi) q[1];
rz(-2.0678287) q[2];
sx q[2];
rz(-1.4364527) q[2];
sx q[2];
rz(1.6430935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.313844) q[1];
sx q[1];
rz(-1.5428468) q[1];
sx q[1];
rz(1.8716783) q[1];
rz(-pi) q[2];
rz(-2.1928093) q[3];
sx q[3];
rz(-1.7025196) q[3];
sx q[3];
rz(-1.458481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.10609145) q[2];
sx q[2];
rz(-1.1320628) q[2];
sx q[2];
rz(2.3512225) q[2];
rz(0.60449374) q[3];
sx q[3];
rz(-2.1166182) q[3];
sx q[3];
rz(0.42009556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6919493) q[0];
sx q[0];
rz(-0.88570166) q[0];
sx q[0];
rz(2.7272136) q[0];
rz(-2.5212133) q[1];
sx q[1];
rz(-1.2216156) q[1];
sx q[1];
rz(1.5827804) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65711439) q[0];
sx q[0];
rz(-0.14095356) q[0];
sx q[0];
rz(-2.4993505) q[0];
x q[1];
rz(-3.0364939) q[2];
sx q[2];
rz(-1.5886279) q[2];
sx q[2];
rz(0.43184973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95629809) q[1];
sx q[1];
rz(-1.5130338) q[1];
sx q[1];
rz(-0.0094780427) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9314194) q[3];
sx q[3];
rz(-2.7322331) q[3];
sx q[3];
rz(1.6334074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15535007) q[2];
sx q[2];
rz(-2.4537931) q[2];
sx q[2];
rz(2.9637994) q[2];
rz(2.6627461) q[3];
sx q[3];
rz(-1.1136473) q[3];
sx q[3];
rz(3.0732885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9947522) q[0];
sx q[0];
rz(-2.1629592) q[0];
sx q[0];
rz(-0.11257182) q[0];
rz(2.5162137) q[1];
sx q[1];
rz(-2.0140779) q[1];
sx q[1];
rz(0.88923997) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93619868) q[0];
sx q[0];
rz(-1.6136716) q[0];
sx q[0];
rz(-1.491719) q[0];
x q[1];
rz(2.8761112) q[2];
sx q[2];
rz(-1.6471787) q[2];
sx q[2];
rz(0.35778174) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4496838) q[1];
sx q[1];
rz(-2.0912716) q[1];
sx q[1];
rz(2.8694105) q[1];
rz(-pi) q[2];
rz(-1.806385) q[3];
sx q[3];
rz(-1.5054678) q[3];
sx q[3];
rz(2.3193086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.858295) q[2];
sx q[2];
rz(-0.12792835) q[2];
sx q[2];
rz(-0.27321401) q[2];
rz(2.9122399) q[3];
sx q[3];
rz(-1.5867686) q[3];
sx q[3];
rz(1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44723085) q[0];
sx q[0];
rz(-2.2706967) q[0];
sx q[0];
rz(-2.2253775) q[0];
rz(-2.9592196) q[1];
sx q[1];
rz(-2.9353751) q[1];
sx q[1];
rz(-1.1590385) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3834476) q[0];
sx q[0];
rz(-0.30320692) q[0];
sx q[0];
rz(2.4571553) q[0];
x q[1];
rz(-2.6508337) q[2];
sx q[2];
rz(-0.10819866) q[2];
sx q[2];
rz(-2.8645017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68284833) q[1];
sx q[1];
rz(-0.35431752) q[1];
sx q[1];
rz(2.0451727) q[1];
rz(-pi) q[2];
rz(2.2044936) q[3];
sx q[3];
rz(-2.1917412) q[3];
sx q[3];
rz(1.6976732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75659043) q[2];
sx q[2];
rz(-0.80005163) q[2];
sx q[2];
rz(-0.315256) q[2];
rz(-2.5473525) q[3];
sx q[3];
rz(-0.88163328) q[3];
sx q[3];
rz(-0.63410223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88157982) q[0];
sx q[0];
rz(-2.4729112) q[0];
sx q[0];
rz(2.9823533) q[0];
rz(1.0881933) q[1];
sx q[1];
rz(-2.2290778) q[1];
sx q[1];
rz(0.27096567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2361901) q[0];
sx q[0];
rz(-1.3199249) q[0];
sx q[0];
rz(-2.0467) q[0];
x q[1];
rz(2.4687211) q[2];
sx q[2];
rz(-1.4904163) q[2];
sx q[2];
rz(-0.56416914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9393255) q[1];
sx q[1];
rz(-1.2346054) q[1];
sx q[1];
rz(2.4656117) q[1];
rz(0.58862092) q[3];
sx q[3];
rz(-0.56097066) q[3];
sx q[3];
rz(-2.4952793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8012041) q[2];
sx q[2];
rz(-0.77216721) q[2];
sx q[2];
rz(2.8188952) q[2];
rz(-3.0835551) q[3];
sx q[3];
rz(-2.3400584) q[3];
sx q[3];
rz(0.29840741) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4651481) q[0];
sx q[0];
rz(-1.5100751) q[0];
sx q[0];
rz(-0.81612192) q[0];
rz(2.8836518) q[1];
sx q[1];
rz(-1.1358658) q[1];
sx q[1];
rz(1.6075016) q[1];
rz(0.77946812) q[2];
sx q[2];
rz(-0.39579724) q[2];
sx q[2];
rz(-0.76553065) q[2];
rz(-2.0097575) q[3];
sx q[3];
rz(-0.96228941) q[3];
sx q[3];
rz(-1.8214772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
