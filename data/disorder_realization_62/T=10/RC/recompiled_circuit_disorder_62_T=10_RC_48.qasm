OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7952607) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(0.46121009) q[0];
rz(1.2167591) q[2];
sx q[2];
rz(-2.3166222) q[2];
sx q[2];
rz(-2.0113457) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8780958) q[1];
sx q[1];
rz(-2.7090008) q[1];
sx q[1];
rz(-2.2357975) q[1];
rz(1.8144242) q[3];
sx q[3];
rz(-1.676179) q[3];
sx q[3];
rz(-1.5308612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(-0.82582981) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.9155496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81235028) q[0];
sx q[0];
rz(-1.7698405) q[0];
sx q[0];
rz(-0.45254405) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3833582) q[2];
sx q[2];
rz(-2.1633534) q[2];
sx q[2];
rz(-0.50298467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1184428) q[1];
sx q[1];
rz(-2.9217302) q[1];
sx q[1];
rz(1.5688194) q[1];
x q[2];
rz(2.1434104) q[3];
sx q[3];
rz(-1.5787573) q[3];
sx q[3];
rz(1.5919459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.9817292) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.8925517) q[0];
rz(0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(2.1121315) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3345966) q[0];
sx q[0];
rz(-0.96141978) q[0];
sx q[0];
rz(0.2683123) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7715893) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(1.5918658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0440812) q[1];
sx q[1];
rz(-1.71002) q[1];
sx q[1];
rz(-0.15686762) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13266487) q[3];
sx q[3];
rz(-2.3972315) q[3];
sx q[3];
rz(-2.6877407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.007894667) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(0.50764817) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(-1.5959651) q[0];
rz(-1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0632616) q[0];
sx q[0];
rz(-0.69561361) q[0];
sx q[0];
rz(2.9177833) q[0];
x q[1];
rz(0.085725351) q[2];
sx q[2];
rz(-2.1282196) q[2];
sx q[2];
rz(0.9312219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1277395) q[1];
sx q[1];
rz(-1.8426367) q[1];
sx q[1];
rz(2.4747162) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32228542) q[3];
sx q[3];
rz(-0.61338131) q[3];
sx q[3];
rz(-1.3410638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76379124) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(0.14136782) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871053) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-0.89486665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7475815) q[0];
sx q[0];
rz(-0.78254875) q[0];
sx q[0];
rz(-2.5400019) q[0];
rz(-1.0084444) q[2];
sx q[2];
rz(-2.214606) q[2];
sx q[2];
rz(1.84562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27069651) q[1];
sx q[1];
rz(-2.0082698) q[1];
sx q[1];
rz(-0.84769627) q[1];
rz(-0.026167913) q[3];
sx q[3];
rz(-0.28655616) q[3];
sx q[3];
rz(-2.9122796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557945) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.3775795) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-0.46498743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12596345) q[0];
sx q[0];
rz(-1.1383801) q[0];
sx q[0];
rz(-2.7727491) q[0];
rz(1.6443406) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(-1.3982915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8556559) q[1];
sx q[1];
rz(-0.93614139) q[1];
sx q[1];
rz(0.059307701) q[1];
x q[2];
rz(1.7808077) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(0.39892808) q[3];
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
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(0.26708189) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12748195) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0359081) q[0];
sx q[0];
rz(-1.3996291) q[0];
sx q[0];
rz(2.8842584) q[0];
rz(1.1440802) q[2];
sx q[2];
rz(-1.0973755) q[2];
sx q[2];
rz(1.13525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1289039) q[1];
sx q[1];
rz(-1.3959937) q[1];
sx q[1];
rz(3.0867982) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3994) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.069313958) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496562) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-0.33777133) q[0];
rz(2.0514964) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(-0.24857323) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0524806) q[0];
sx q[0];
rz(-0.57894527) q[0];
sx q[0];
rz(-0.46778932) q[0];
rz(-pi) q[1];
rz(0.82069223) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(1.5936268) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0174745) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(0.89213051) q[1];
rz(-pi) q[2];
rz(-1.9931273) q[3];
sx q[3];
rz(-0.18581192) q[3];
sx q[3];
rz(-1.2687792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(0.48197204) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(1.1530676) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865737) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(-1.823002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5134207) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(-2.0128987) q[0];
x q[1];
rz(-2.4039688) q[2];
sx q[2];
rz(-2.3361623) q[2];
sx q[2];
rz(-1.1951624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2576037) q[1];
sx q[1];
rz(-0.89234967) q[1];
sx q[1];
rz(0.18737327) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3913279) q[3];
sx q[3];
rz(-1.8599469) q[3];
sx q[3];
rz(0.17444785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(2.7424157) q[2];
rz(0.88360751) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(1.261196) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8226782) q[0];
sx q[0];
rz(-0.42352522) q[0];
sx q[0];
rz(0.25869297) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.10896) q[2];
sx q[2];
rz(-2.543078) q[2];
sx q[2];
rz(-1.061071) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1120158) q[1];
sx q[1];
rz(-0.50260168) q[1];
sx q[1];
rz(0.13346787) q[1];
rz(-pi) q[2];
rz(2.946978) q[3];
sx q[3];
rz(-0.6932887) q[3];
sx q[3];
rz(0.64893901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(2.5718001) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-0.56263721) q[3];
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
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(-0.60824153) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(1.6115887) q[2];
sx q[2];
rz(-0.57208021) q[2];
sx q[2];
rz(-0.013442599) q[2];
rz(1.431987) q[3];
sx q[3];
rz(-2.1574253) q[3];
sx q[3];
rz(0.16711259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
