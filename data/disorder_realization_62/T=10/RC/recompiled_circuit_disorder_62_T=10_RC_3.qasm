OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(-2.9641889) q[0];
sx q[0];
rz(2.0071964) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21270574) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(1.5972208) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9248336) q[2];
sx q[2];
rz(-0.82497043) q[2];
sx q[2];
rz(1.1302469) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2156148) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(-1.2222626) q[1];
x q[2];
rz(1.1575559) q[3];
sx q[3];
rz(-2.8765656) q[3];
sx q[3];
rz(-2.7812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87876451) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-0.051068548) q[2];
rz(2.5845394) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(-2.5449975) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(1.9155496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37218371) q[0];
sx q[0];
rz(-2.6499977) q[0];
sx q[0];
rz(2.7093637) q[0];
rz(-2.3833582) q[2];
sx q[2];
rz(-0.97823921) q[2];
sx q[2];
rz(2.638608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.021124161) q[1];
sx q[1];
rz(-1.3509343) q[1];
sx q[1];
rz(3.1411509) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5854884) q[3];
sx q[3];
rz(-0.5726632) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(2.8105695) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.249041) q[0];
rz(0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(2.1121315) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345966) q[0];
sx q[0];
rz(-2.1801729) q[0];
sx q[0];
rz(-2.8732804) q[0];
x q[1];
rz(-1.2842032) q[2];
sx q[2];
rz(-2.2006052) q[2];
sx q[2];
rz(-2.0558002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.894695) q[1];
sx q[1];
rz(-0.20935911) q[1];
sx q[1];
rz(-0.73114242) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13266487) q[3];
sx q[3];
rz(-0.74436114) q[3];
sx q[3];
rz(0.45385195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(-2.6339445) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.1631789) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65748173) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(-1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31947485) q[0];
sx q[0];
rz(-1.7135156) q[0];
sx q[0];
rz(-0.6832173) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0117202) q[2];
sx q[2];
rz(-1.6435197) q[2];
sx q[2];
rz(0.59414547) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1277395) q[1];
sx q[1];
rz(-1.8426367) q[1];
sx q[1];
rz(0.66687648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3514148) q[3];
sx q[3];
rz(-0.99321584) q[3];
sx q[3];
rz(-2.1882309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871053) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8653523) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(0.68666896) q[0];
x q[1];
rz(2.4159555) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(2.5051136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85256165) q[1];
sx q[1];
rz(-2.3173996) q[1];
sx q[1];
rz(-2.186071) q[1];
rz(-pi) q[2];
rz(-1.5785061) q[3];
sx q[3];
rz(-1.8572516) q[3];
sx q[3];
rz(0.25659284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.3775795) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-2.6766052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2705921) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(0.90765783) q[0];
rz(-pi) q[1];
rz(-1.6443406) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(1.3982915) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8556559) q[1];
sx q[1];
rz(-2.2054513) q[1];
sx q[1];
rz(3.082285) q[1];
x q[2];
rz(1.3607849) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(2.7426646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0141107) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(-1.4200462) q[0];
rz(-3.1177915) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0359081) q[0];
sx q[0];
rz(-1.7419635) q[0];
sx q[0];
rz(-0.25733421) q[0];
rz(-pi) q[1];
rz(-2.4620352) q[2];
sx q[2];
rz(-2.5153) q[2];
sx q[2];
rz(-1.9192413) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29282001) q[1];
sx q[1];
rz(-0.18310586) q[1];
sx q[1];
rz(1.2700901) q[1];
rz(-pi) q[2];
rz(2.3994) q[3];
sx q[3];
rz(-0.67897292) q[3];
sx q[3];
rz(1.2364482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-2.0689266) q[2];
rz(0.32564751) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(1.1919365) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(2.8038213) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(2.8930194) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2598495) q[0];
sx q[0];
rz(-1.3215085) q[0];
sx q[0];
rz(-0.52815204) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48034251) q[2];
sx q[2];
rz(-1.1104274) q[2];
sx q[2];
rz(0.72887052) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92330248) q[1];
sx q[1];
rz(-2.3308838) q[1];
sx q[1];
rz(0.87358012) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9931273) q[3];
sx q[3];
rz(-2.9557807) q[3];
sx q[3];
rz(1.8728135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(-1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(0.28276643) q[0];
rz(-2.4400318) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.3185906) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2254927) q[0];
sx q[0];
rz(-0.44888228) q[0];
sx q[0];
rz(-1.7569957) q[0];
rz(2.4039688) q[2];
sx q[2];
rz(-0.80543033) q[2];
sx q[2];
rz(-1.1951624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1773771) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(1.3436951) q[1];
rz(-pi) q[2];
rz(0.75026476) q[3];
sx q[3];
rz(-1.8599469) q[3];
sx q[3];
rz(2.9671448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(2.7424157) q[2];
rz(-0.88360751) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.2214899) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5426853) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.8803966) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529918) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(-2.7306042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59826675) q[2];
sx q[2];
rz(-1.5524128) q[2];
sx q[2];
rz(-0.53668864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8775455) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(-1.4977786) q[1];
x q[2];
rz(-0.68393771) q[3];
sx q[3];
rz(-1.6947019) q[3];
sx q[3];
rz(-2.370196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8284843) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-0.026253168) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(1.7096056) q[3];
sx q[3];
rz(-0.98416735) q[3];
sx q[3];
rz(-2.9744801) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
