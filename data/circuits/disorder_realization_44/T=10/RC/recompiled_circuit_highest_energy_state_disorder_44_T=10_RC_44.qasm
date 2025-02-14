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
rz(1.1062082) q[0];
sx q[0];
rz(-2.455403) q[0];
sx q[0];
rz(-1.9880779) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(0.44841132) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2758688) q[0];
sx q[0];
rz(-0.65058904) q[0];
sx q[0];
rz(-0.087819789) q[0];
rz(-pi) q[1];
rz(-1.6719867) q[2];
sx q[2];
rz(-1.8253583) q[2];
sx q[2];
rz(0.86092608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37473512) q[1];
sx q[1];
rz(-2.285706) q[1];
sx q[1];
rz(-0.84918569) q[1];
rz(-0.5011933) q[3];
sx q[3];
rz(-1.5287011) q[3];
sx q[3];
rz(0.30748366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6260234) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-2.5173729) q[2];
rz(2.7624687) q[3];
sx q[3];
rz(-2.1513394) q[3];
sx q[3];
rz(2.8930801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(-0.17495951) q[0];
sx q[0];
rz(-2.3269854) q[0];
sx q[0];
rz(2.1061184) q[0];
rz(-1.8677208) q[1];
sx q[1];
rz(-2.404411) q[1];
sx q[1];
rz(0.48286352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55601701) q[0];
sx q[0];
rz(-3.034798) q[0];
sx q[0];
rz(0.44321816) q[0];
rz(-0.61699485) q[2];
sx q[2];
rz(-0.20295396) q[2];
sx q[2];
rz(-2.647612) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4299996) q[1];
sx q[1];
rz(-1.5614127) q[1];
sx q[1];
rz(2.0026155) q[1];
rz(-2.9008174) q[3];
sx q[3];
rz(-2.8117883) q[3];
sx q[3];
rz(2.2036116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1172993) q[2];
sx q[2];
rz(-1.6397986) q[2];
sx q[2];
rz(-1.2705605) q[2];
rz(0.67000669) q[3];
sx q[3];
rz(-2.5195401) q[3];
sx q[3];
rz(0.89699927) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73827493) q[0];
sx q[0];
rz(-2.6955695) q[0];
sx q[0];
rz(0.73905149) q[0];
rz(-2.1801379) q[1];
sx q[1];
rz(-0.32197222) q[1];
sx q[1];
rz(0.75698537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3630545) q[0];
sx q[0];
rz(-0.86150384) q[0];
sx q[0];
rz(-0.34687931) q[0];
x q[1];
rz(0.69125533) q[2];
sx q[2];
rz(-1.2757841) q[2];
sx q[2];
rz(-1.472109) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0464335) q[1];
sx q[1];
rz(-0.94645247) q[1];
sx q[1];
rz(-1.4051564) q[1];
x q[2];
rz(1.3540165) q[3];
sx q[3];
rz(-0.31459537) q[3];
sx q[3];
rz(1.6198688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5028533) q[2];
sx q[2];
rz(-0.91705489) q[2];
sx q[2];
rz(2.4364831) q[2];
rz(-0.32143587) q[3];
sx q[3];
rz(-1.0175846) q[3];
sx q[3];
rz(0.41079918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.131677) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(1.2257082) q[0];
rz(1.2340087) q[1];
sx q[1];
rz(-1.0154513) q[1];
sx q[1];
rz(-0.066224901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4915337) q[0];
sx q[0];
rz(-0.26187944) q[0];
sx q[0];
rz(-1.2175548) q[0];
rz(-pi) q[1];
x q[1];
rz(2.340256) q[2];
sx q[2];
rz(-1.7883915) q[2];
sx q[2];
rz(-0.8720397) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.855763) q[1];
sx q[1];
rz(-2.0501839) q[1];
sx q[1];
rz(-2.781032) q[1];
rz(-pi) q[2];
rz(-0.31853557) q[3];
sx q[3];
rz(-0.38357601) q[3];
sx q[3];
rz(0.84795241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0951198) q[2];
sx q[2];
rz(-2.1472223) q[2];
sx q[2];
rz(2.7009916) q[2];
rz(-2.1353841) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(2.3073176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7202268) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(1.45586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86690534) q[0];
sx q[0];
rz(-0.84466776) q[0];
sx q[0];
rz(-2.9466036) q[0];
rz(-pi) q[1];
rz(-2.9142889) q[2];
sx q[2];
rz(-1.6844201) q[2];
sx q[2];
rz(-2.794968) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2537186) q[1];
sx q[1];
rz(-1.5504596) q[1];
sx q[1];
rz(0.65171297) q[1];
x q[2];
rz(2.9221228) q[3];
sx q[3];
rz(-0.61958909) q[3];
sx q[3];
rz(-2.0989024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8684034) q[2];
sx q[2];
rz(-1.5634147) q[2];
sx q[2];
rz(-0.91147649) q[2];
rz(-2.2677926) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(-3.1349283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28354302) q[0];
sx q[0];
rz(-0.25074211) q[0];
sx q[0];
rz(0.13033303) q[0];
rz(0.48769543) q[1];
sx q[1];
rz(-2.6002488) q[1];
sx q[1];
rz(0.74388751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.795563) q[0];
sx q[0];
rz(-1.6206121) q[0];
sx q[0];
rz(-1.6532142) q[0];
rz(0.13915542) q[2];
sx q[2];
rz(-0.94991131) q[2];
sx q[2];
rz(-0.54229743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6106657) q[1];
sx q[1];
rz(-2.1082889) q[1];
sx q[1];
rz(-2.8100138) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3906995) q[3];
sx q[3];
rz(-1.0922179) q[3];
sx q[3];
rz(-3.0949926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1193715) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(0.96088299) q[2];
rz(-0.39572257) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(-3.0254288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845067) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(1.047026) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.2774757) q[1];
sx q[1];
rz(0.39271694) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1395124) q[0];
sx q[0];
rz(-0.89186984) q[0];
sx q[0];
rz(0.2261136) q[0];
rz(1.4899859) q[2];
sx q[2];
rz(-2.5210292) q[2];
sx q[2];
rz(-2.0106237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9153292) q[1];
sx q[1];
rz(-1.8460423) q[1];
sx q[1];
rz(-1.5987426) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3052528) q[3];
sx q[3];
rz(-1.8887541) q[3];
sx q[3];
rz(-1.6443271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.087223209) q[2];
sx q[2];
rz(-1.1944218) q[2];
sx q[2];
rz(-1.3090022) q[2];
rz(1.7878112) q[3];
sx q[3];
rz(-1.9446707) q[3];
sx q[3];
rz(0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29276174) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(2.8344179) q[0];
rz(-1.81709) q[1];
sx q[1];
rz(-0.38506404) q[1];
sx q[1];
rz(-0.90075341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9884315) q[0];
sx q[0];
rz(-1.5539196) q[0];
sx q[0];
rz(0.26217006) q[0];
rz(-pi) q[1];
rz(2.1395166) q[2];
sx q[2];
rz(-0.7373215) q[2];
sx q[2];
rz(-1.3739283) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72462725) q[1];
sx q[1];
rz(-2.6748383) q[1];
sx q[1];
rz(0.80893597) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19260223) q[3];
sx q[3];
rz(-0.64036548) q[3];
sx q[3];
rz(-2.7964724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8073392) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(2.5267498) q[2];
rz(1.0963415) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(-2.443327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7480302) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(2.6322741) q[0];
rz(2.9604984) q[1];
sx q[1];
rz(-1.4053922) q[1];
sx q[1];
rz(-0.96955713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90617243) q[0];
sx q[0];
rz(-2.3114822) q[0];
sx q[0];
rz(-1.7801911) q[0];
x q[1];
rz(-2.0779586) q[2];
sx q[2];
rz(-1.5263057) q[2];
sx q[2];
rz(-2.8756093) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4050715) q[1];
sx q[1];
rz(-0.35686359) q[1];
sx q[1];
rz(0.31330152) q[1];
rz(-0.35004444) q[3];
sx q[3];
rz(-2.2681103) q[3];
sx q[3];
rz(0.75087912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2481044) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(1.0815557) q[2];
rz(-0.23165101) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(-2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.76081) q[0];
sx q[0];
rz(-0.2247227) q[0];
sx q[0];
rz(2.494452) q[0];
rz(0.45516792) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(-0.761935) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021558048) q[0];
sx q[0];
rz(-1.8592478) q[0];
sx q[0];
rz(2.4214273) q[0];
x q[1];
rz(-2.5013431) q[2];
sx q[2];
rz(-0.92776042) q[2];
sx q[2];
rz(-1.4852448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0207723) q[1];
sx q[1];
rz(-1.5365531) q[1];
sx q[1];
rz(-2.4832151) q[1];
rz(-pi) q[2];
rz(1.5555218) q[3];
sx q[3];
rz(-0.46746436) q[3];
sx q[3];
rz(-1.8031507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.053293856) q[2];
sx q[2];
rz(-0.21685728) q[2];
sx q[2];
rz(-2.2376412) q[2];
rz(1.5229185) q[3];
sx q[3];
rz(-1.0205597) q[3];
sx q[3];
rz(-1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91697964) q[0];
sx q[0];
rz(-1.1876748) q[0];
sx q[0];
rz(-1.351958) q[0];
rz(3.0147973) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(2.0595111) q[2];
sx q[2];
rz(-0.50177028) q[2];
sx q[2];
rz(-0.72015464) q[2];
rz(-0.0062777304) q[3];
sx q[3];
rz(-2.7115887) q[3];
sx q[3];
rz(-0.030621519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
