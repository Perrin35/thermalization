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
rz(0.25207818) q[0];
sx q[0];
rz(3.7171465) q[0];
sx q[0];
rz(9.3024749) q[0];
rz(1.1072493) q[1];
sx q[1];
rz(-0.9613494) q[1];
sx q[1];
rz(-0.038318757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81403804) q[0];
sx q[0];
rz(-0.71056847) q[0];
sx q[0];
rz(2.1585805) q[0];
rz(-pi) q[1];
rz(-3.0228268) q[2];
sx q[2];
rz(-1.6890172) q[2];
sx q[2];
rz(1.880065) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1120834) q[1];
sx q[1];
rz(-1.2473628) q[1];
sx q[1];
rz(2.1548163) q[1];
x q[2];
rz(-1.5502717) q[3];
sx q[3];
rz(-1.3234659) q[3];
sx q[3];
rz(1.0350682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4965839) q[2];
sx q[2];
rz(-1.714548) q[2];
sx q[2];
rz(-1.6928147) q[2];
rz(2.951176) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(-0.14934389) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377624) q[0];
sx q[0];
rz(-1.2685403) q[0];
sx q[0];
rz(2.8835836) q[0];
rz(1.5915271) q[1];
sx q[1];
rz(-1.9015046) q[1];
sx q[1];
rz(2.9749427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7014871) q[0];
sx q[0];
rz(-0.99033725) q[0];
sx q[0];
rz(1.6850182) q[0];
rz(-pi) q[1];
rz(-2.7562253) q[2];
sx q[2];
rz(-1.8358942) q[2];
sx q[2];
rz(-0.4601269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38063654) q[1];
sx q[1];
rz(-1.6957307) q[1];
sx q[1];
rz(-2.4255468) q[1];
rz(-0.54661481) q[3];
sx q[3];
rz(-1.0738465) q[3];
sx q[3];
rz(0.95595804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0717281) q[2];
sx q[2];
rz(-2.0191777) q[2];
sx q[2];
rz(1.9094763) q[2];
rz(2.2169436) q[3];
sx q[3];
rz(-0.93622127) q[3];
sx q[3];
rz(2.0360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.449618) q[0];
sx q[0];
rz(-1.6716577) q[0];
sx q[0];
rz(0.0019419226) q[0];
rz(3.107403) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(1.5431822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6337834) q[0];
sx q[0];
rz(-2.6144321) q[0];
sx q[0];
rz(-2.6109004) q[0];
rz(-pi) q[1];
rz(2.4056881) q[2];
sx q[2];
rz(-1.8712448) q[2];
sx q[2];
rz(1.7395333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.10101) q[1];
sx q[1];
rz(-1.2173182) q[1];
sx q[1];
rz(1.6124875) q[1];
x q[2];
rz(-2.8974124) q[3];
sx q[3];
rz(-2.2043318) q[3];
sx q[3];
rz(1.4118495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89692846) q[2];
sx q[2];
rz(-0.43411532) q[2];
sx q[2];
rz(1.0373235) q[2];
rz(-1.1050998) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(1.1387811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8741375) q[0];
sx q[0];
rz(-2.656811) q[0];
sx q[0];
rz(-1.7304035) q[0];
rz(2.5022068) q[1];
sx q[1];
rz(-1.4566028) q[1];
sx q[1];
rz(0.04714084) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9586724) q[0];
sx q[0];
rz(-2.483568) q[0];
sx q[0];
rz(0.47112314) q[0];
rz(-pi) q[1];
rz(2.7732255) q[2];
sx q[2];
rz(-0.50058156) q[2];
sx q[2];
rz(2.6427302) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5410536) q[1];
sx q[1];
rz(-2.4324634) q[1];
sx q[1];
rz(1.0102655) q[1];
rz(-pi) q[2];
rz(-1.9023667) q[3];
sx q[3];
rz(-0.75115381) q[3];
sx q[3];
rz(1.0855261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3186657) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(1.218943) q[2];
rz(-2.8259891) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3217992) q[0];
sx q[0];
rz(-2.5892374) q[0];
sx q[0];
rz(-0.34859443) q[0];
rz(-1.2527342) q[1];
sx q[1];
rz(-1.9786973) q[1];
sx q[1];
rz(1.3016275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5900629) q[0];
sx q[0];
rz(-0.7669391) q[0];
sx q[0];
rz(1.4566896) q[0];
rz(-pi) q[1];
rz(2.581004) q[2];
sx q[2];
rz(-1.9288256) q[2];
sx q[2];
rz(2.5303417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5432413) q[1];
sx q[1];
rz(-1.390402) q[1];
sx q[1];
rz(1.3460657) q[1];
x q[2];
rz(0.35195203) q[3];
sx q[3];
rz(-1.0478596) q[3];
sx q[3];
rz(0.74229303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9584413) q[2];
sx q[2];
rz(-0.30854598) q[2];
sx q[2];
rz(1.4754254) q[2];
rz(-0.685855) q[3];
sx q[3];
rz(-2.1619022) q[3];
sx q[3];
rz(0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0297861) q[0];
sx q[0];
rz(-2.989558) q[0];
sx q[0];
rz(1.6492122) q[0];
rz(-2.1624508) q[1];
sx q[1];
rz(-1.6056332) q[1];
sx q[1];
rz(1.917256) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79681841) q[0];
sx q[0];
rz(-2.9412656) q[0];
sx q[0];
rz(-2.9873559) q[0];
rz(-pi) q[1];
rz(-0.54927214) q[2];
sx q[2];
rz(-1.7584929) q[2];
sx q[2];
rz(0.21085462) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.853272) q[1];
sx q[1];
rz(-1.558312) q[1];
sx q[1];
rz(-2.1846143) q[1];
rz(-0.94131366) q[3];
sx q[3];
rz(-1.8555897) q[3];
sx q[3];
rz(2.3465057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7225723) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(-2.2255619) q[2];
rz(-1.3845059) q[3];
sx q[3];
rz(-1.2576831) q[3];
sx q[3];
rz(-2.4345583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0093805669) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(0.95426553) q[0];
rz(-2.7047899) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(-0.11016914) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018579114) q[0];
sx q[0];
rz(-1.8948104) q[0];
sx q[0];
rz(-1.8229681) q[0];
rz(0.61715257) q[2];
sx q[2];
rz(-0.79884702) q[2];
sx q[2];
rz(-1.3399194) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72790775) q[1];
sx q[1];
rz(-1.2251108) q[1];
sx q[1];
rz(3.137008) q[1];
rz(-pi) q[2];
rz(0.25400587) q[3];
sx q[3];
rz(-1.0830823) q[3];
sx q[3];
rz(0.28764492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7097912) q[2];
sx q[2];
rz(-3.1352037) q[2];
sx q[2];
rz(-2.1978281) q[2];
rz(-2.5778095) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(1.110466) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540045) q[0];
sx q[0];
rz(-2.8546951) q[0];
sx q[0];
rz(0.34550825) q[0];
rz(-1.0954789) q[1];
sx q[1];
rz(-1.6637207) q[1];
sx q[1];
rz(1.5617721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0559383) q[0];
sx q[0];
rz(-1.5393747) q[0];
sx q[0];
rz(2.3532193) q[0];
x q[1];
rz(1.4905632) q[2];
sx q[2];
rz(-0.65872619) q[2];
sx q[2];
rz(-1.4874489) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53198481) q[1];
sx q[1];
rz(-0.72276211) q[1];
sx q[1];
rz(-1.0015798) q[1];
x q[2];
rz(-1.5388266) q[3];
sx q[3];
rz(-1.073146) q[3];
sx q[3];
rz(-0.15807334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3788508) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(1.4307107) q[2];
rz(2.6070969) q[3];
sx q[3];
rz(-1.675324) q[3];
sx q[3];
rz(2.1889595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6119824) q[0];
sx q[0];
rz(-0.069267608) q[0];
sx q[0];
rz(1.8310504) q[0];
rz(-0.88690859) q[1];
sx q[1];
rz(-0.84019089) q[1];
sx q[1];
rz(-1.8720253) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55063215) q[0];
sx q[0];
rz(-2.7661588) q[0];
sx q[0];
rz(-1.9272789) q[0];
x q[1];
rz(1.1940895) q[2];
sx q[2];
rz(-2.7393638) q[2];
sx q[2];
rz(-2.5144983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15358217) q[1];
sx q[1];
rz(-1.8021823) q[1];
sx q[1];
rz(-1.2910976) q[1];
x q[2];
rz(0.2036152) q[3];
sx q[3];
rz(-1.5316267) q[3];
sx q[3];
rz(2.7050582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5294007) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(-0.23294918) q[2];
rz(1.285078) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(-0.3860093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9832298) q[0];
sx q[0];
rz(-2.5322999) q[0];
sx q[0];
rz(3.0443211) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(2.9248617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6410206) q[0];
sx q[0];
rz(-1.4002698) q[0];
sx q[0];
rz(3.0855623) q[0];
rz(-1.6172269) q[2];
sx q[2];
rz(-1.7123607) q[2];
sx q[2];
rz(2.9449449) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3236802) q[1];
sx q[1];
rz(-0.92880946) q[1];
sx q[1];
rz(0.85425185) q[1];
rz(2.7085767) q[3];
sx q[3];
rz(-2.4320514) q[3];
sx q[3];
rz(2.7482928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8486166) q[2];
sx q[2];
rz(-2.8605707) q[2];
sx q[2];
rz(2.6463553) q[2];
rz(-1.5826591) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(-0.16286477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0804629) q[0];
sx q[0];
rz(-1.6292138) q[0];
sx q[0];
rz(-0.52393352) q[0];
rz(-1.0864661) q[1];
sx q[1];
rz(-2.6230984) q[1];
sx q[1];
rz(-2.7883504) q[1];
rz(-1.6025966) q[2];
sx q[2];
rz(-1.5780137) q[2];
sx q[2];
rz(0.93102166) q[2];
rz(-0.85862715) q[3];
sx q[3];
rz(-1.5024019) q[3];
sx q[3];
rz(2.0988322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
