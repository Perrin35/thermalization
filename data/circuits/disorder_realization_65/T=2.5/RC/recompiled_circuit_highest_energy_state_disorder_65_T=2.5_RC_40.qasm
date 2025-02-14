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
rz(0.12198099) q[0];
sx q[0];
rz(3.0966336) q[0];
sx q[0];
rz(9.3293204) q[0];
rz(1.8010315) q[1];
sx q[1];
rz(3.0756693) q[1];
sx q[1];
rz(8.8311721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5779764) q[0];
sx q[0];
rz(-1.9270183) q[0];
sx q[0];
rz(2.6249159) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6907211) q[2];
sx q[2];
rz(-1.7210135) q[2];
sx q[2];
rz(0.37301979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.846147) q[1];
sx q[1];
rz(-2.4762003) q[1];
sx q[1];
rz(-0.34122463) q[1];
rz(2.8659561) q[3];
sx q[3];
rz(-2.079981) q[3];
sx q[3];
rz(1.0643626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.060071271) q[2];
sx q[2];
rz(-1.8895443) q[2];
sx q[2];
rz(-2.579465) q[2];
rz(-2.8020322) q[3];
sx q[3];
rz(-1.5226676) q[3];
sx q[3];
rz(-1.08574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.958441) q[0];
sx q[0];
rz(-0.14587942) q[0];
sx q[0];
rz(-1.1613783) q[0];
rz(0.24102744) q[1];
sx q[1];
rz(-0.82254219) q[1];
sx q[1];
rz(2.8384812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2456692) q[0];
sx q[0];
rz(-1.5035986) q[0];
sx q[0];
rz(1.3724061) q[0];
x q[1];
rz(-2.9955818) q[2];
sx q[2];
rz(-1.6939298) q[2];
sx q[2];
rz(1.6347803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0063340291) q[1];
sx q[1];
rz(-2.5016682) q[1];
sx q[1];
rz(-2.695822) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47789033) q[3];
sx q[3];
rz(-1.9955561) q[3];
sx q[3];
rz(0.44179824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8384398) q[2];
sx q[2];
rz(-1.2829605) q[2];
sx q[2];
rz(-2.8008833) q[2];
rz(-0.43863145) q[3];
sx q[3];
rz(-0.80629587) q[3];
sx q[3];
rz(0.058569245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230147) q[0];
sx q[0];
rz(-0.61841643) q[0];
sx q[0];
rz(-0.83786905) q[0];
rz(-2.2184929) q[1];
sx q[1];
rz(-1.9337312) q[1];
sx q[1];
rz(2.7941678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9073931) q[0];
sx q[0];
rz(-0.1009909) q[0];
sx q[0];
rz(-1.9491461) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4342257) q[2];
sx q[2];
rz(-2.8828808) q[2];
sx q[2];
rz(1.9345891) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8953462) q[1];
sx q[1];
rz(-2.0450174) q[1];
sx q[1];
rz(-0.62690027) q[1];
x q[2];
rz(2.4597563) q[3];
sx q[3];
rz(-0.41636514) q[3];
sx q[3];
rz(-1.8133461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0190987) q[2];
sx q[2];
rz(-0.118003) q[2];
sx q[2];
rz(0.72354358) q[2];
rz(1.8540234) q[3];
sx q[3];
rz(-2.5836594) q[3];
sx q[3];
rz(3.0961228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36640722) q[0];
sx q[0];
rz(-0.43908304) q[0];
sx q[0];
rz(-1.084569) q[0];
rz(-1.1868125) q[1];
sx q[1];
rz(-1.908952) q[1];
sx q[1];
rz(-1.4694227) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050768269) q[0];
sx q[0];
rz(-1.5406797) q[0];
sx q[0];
rz(-1.5876905) q[0];
rz(1.9590366) q[2];
sx q[2];
rz(-1.2080384) q[2];
sx q[2];
rz(1.658939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9778825) q[1];
sx q[1];
rz(-0.61874604) q[1];
sx q[1];
rz(2.2676022) q[1];
rz(-pi) q[2];
rz(0.74923781) q[3];
sx q[3];
rz(-2.6686274) q[3];
sx q[3];
rz(1.7974896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0161418) q[2];
sx q[2];
rz(-0.90105385) q[2];
sx q[2];
rz(2.2060642) q[2];
rz(2.9901796) q[3];
sx q[3];
rz(-0.97984034) q[3];
sx q[3];
rz(-2.6305731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147515) q[0];
sx q[0];
rz(-1.5837357) q[0];
sx q[0];
rz(2.3874808) q[0];
rz(-0.61112815) q[1];
sx q[1];
rz(-1.1437623) q[1];
sx q[1];
rz(-0.67620826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0128764) q[0];
sx q[0];
rz(-2.7749847) q[0];
sx q[0];
rz(-1.104273) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36415139) q[2];
sx q[2];
rz(-1.0478579) q[2];
sx q[2];
rz(2.5464818) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0724746) q[1];
sx q[1];
rz(-2.9058911) q[1];
sx q[1];
rz(-0.75171296) q[1];
x q[2];
rz(1.4884454) q[3];
sx q[3];
rz(-2.1550278) q[3];
sx q[3];
rz(1.0109993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0040697441) q[2];
sx q[2];
rz(-0.71228945) q[2];
sx q[2];
rz(0.75312692) q[2];
rz(-0.5641886) q[3];
sx q[3];
rz(-1.1078395) q[3];
sx q[3];
rz(-2.3851725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18993987) q[0];
sx q[0];
rz(-2.6615182) q[0];
sx q[0];
rz(-2.9204364) q[0];
rz(0.13825522) q[1];
sx q[1];
rz(-1.4555376) q[1];
sx q[1];
rz(2.5855605) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5297896) q[0];
sx q[0];
rz(-1.3767813) q[0];
sx q[0];
rz(0.24873269) q[0];
rz(-2.7660335) q[2];
sx q[2];
rz(-1.9347768) q[2];
sx q[2];
rz(3.08422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7307463) q[1];
sx q[1];
rz(-1.4173943) q[1];
sx q[1];
rz(-1.549841) q[1];
rz(-pi) q[2];
rz(-1.3294499) q[3];
sx q[3];
rz(-2.6785319) q[3];
sx q[3];
rz(-1.9928513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9941462) q[2];
sx q[2];
rz(-0.7615971) q[2];
sx q[2];
rz(0.93696326) q[2];
rz(-2.7312036) q[3];
sx q[3];
rz(-2.161721) q[3];
sx q[3];
rz(2.4858937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3963102) q[0];
sx q[0];
rz(-2.2672125) q[0];
sx q[0];
rz(2.1790047) q[0];
rz(0.67086041) q[1];
sx q[1];
rz(-0.52898359) q[1];
sx q[1];
rz(2.3765391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7350525) q[0];
sx q[0];
rz(-2.708754) q[0];
sx q[0];
rz(1.1821683) q[0];
x q[1];
rz(-0.069326055) q[2];
sx q[2];
rz(-1.6717807) q[2];
sx q[2];
rz(-0.19420964) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.194696) q[1];
sx q[1];
rz(-1.0651154) q[1];
sx q[1];
rz(2.6891788) q[1];
x q[2];
rz(1.9387127) q[3];
sx q[3];
rz(-1.4548317) q[3];
sx q[3];
rz(2.396317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.700231) q[2];
sx q[2];
rz(-0.92741489) q[2];
sx q[2];
rz(2.3564763) q[2];
rz(-0.49078068) q[3];
sx q[3];
rz(-1.1622585) q[3];
sx q[3];
rz(-0.47500113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15071507) q[0];
sx q[0];
rz(-0.10902037) q[0];
sx q[0];
rz(0.14608598) q[0];
rz(0.27169216) q[1];
sx q[1];
rz(-0.79333317) q[1];
sx q[1];
rz(-0.21154107) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5141895) q[0];
sx q[0];
rz(-2.3334921) q[0];
sx q[0];
rz(1.3405728) q[0];
rz(-0.11577932) q[2];
sx q[2];
rz(-1.2010788) q[2];
sx q[2];
rz(0.055495128) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.614173) q[1];
sx q[1];
rz(-2.895311) q[1];
sx q[1];
rz(-0.28913943) q[1];
rz(2.170874) q[3];
sx q[3];
rz(-1.479199) q[3];
sx q[3];
rz(-0.84176871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3286256) q[2];
sx q[2];
rz(-0.93572891) q[2];
sx q[2];
rz(2.3198371) q[2];
rz(1.3385319) q[3];
sx q[3];
rz(-2.0526363) q[3];
sx q[3];
rz(0.72430044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61447918) q[0];
sx q[0];
rz(-2.4452657) q[0];
sx q[0];
rz(1.3592199) q[0];
rz(0.95787734) q[1];
sx q[1];
rz(-0.33708894) q[1];
sx q[1];
rz(-0.27761308) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2569297) q[0];
sx q[0];
rz(-2.1799934) q[0];
sx q[0];
rz(-0.67237396) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3657655) q[2];
sx q[2];
rz(-1.3444573) q[2];
sx q[2];
rz(0.090858484) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23944868) q[1];
sx q[1];
rz(-0.43593299) q[1];
sx q[1];
rz(2.8365342) q[1];
x q[2];
rz(1.3986392) q[3];
sx q[3];
rz(-1.1973698) q[3];
sx q[3];
rz(2.7368054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17310384) q[2];
sx q[2];
rz(-2.5471881) q[2];
sx q[2];
rz(-1.8572726) q[2];
rz(-0.77749085) q[3];
sx q[3];
rz(-2.5992664) q[3];
sx q[3];
rz(0.29977453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(3.1182627) q[0];
sx q[0];
rz(-1.6602004) q[0];
sx q[0];
rz(0.0059286038) q[0];
rz(1.3395576) q[1];
sx q[1];
rz(-2.1052994) q[1];
sx q[1];
rz(-2.6609227) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4673432) q[0];
sx q[0];
rz(-0.021585781) q[0];
sx q[0];
rz(0.72294803) q[0];
x q[1];
rz(1.0239915) q[2];
sx q[2];
rz(-0.725774) q[2];
sx q[2];
rz(-0.99601907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0379767) q[1];
sx q[1];
rz(-2.3876752) q[1];
sx q[1];
rz(2.5186954) q[1];
rz(-pi) q[2];
rz(1.8371498) q[3];
sx q[3];
rz(-1.7031809) q[3];
sx q[3];
rz(1.3416106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9504451) q[2];
sx q[2];
rz(-0.55926776) q[2];
sx q[2];
rz(-0.17678235) q[2];
rz(2.2117129) q[3];
sx q[3];
rz(-1.1888489) q[3];
sx q[3];
rz(1.0298347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4491691) q[0];
sx q[0];
rz(-1.7208736) q[0];
sx q[0];
rz(2.0509913) q[0];
rz(-2.069166) q[1];
sx q[1];
rz(-1.4785531) q[1];
sx q[1];
rz(-1.4428152) q[1];
rz(2.6497447) q[2];
sx q[2];
rz(-1.576423) q[2];
sx q[2];
rz(1.5058422) q[2];
rz(-2.5363793) q[3];
sx q[3];
rz(-1.6463668) q[3];
sx q[3];
rz(2.2201408) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
