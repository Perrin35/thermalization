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
rz(0.81646252) q[0];
sx q[0];
rz(3.2433885) q[0];
sx q[0];
rz(9.9643702) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(0.53909477) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70654725) q[0];
sx q[0];
rz(-2.6994978) q[0];
sx q[0];
rz(-0.059414359) q[0];
rz(-1.4194054) q[2];
sx q[2];
rz(-1.9266204) q[2];
sx q[2];
rz(0.21499888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0338933) q[1];
sx q[1];
rz(-1.7105719) q[1];
sx q[1];
rz(1.2453394) q[1];
rz(0.92181262) q[3];
sx q[3];
rz(-1.6361437) q[3];
sx q[3];
rz(1.9416888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3794136) q[2];
sx q[2];
rz(-0.87612408) q[2];
sx q[2];
rz(-1.1890821) q[2];
rz(-1.1842229) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9963843) q[0];
sx q[0];
rz(-1.5518016) q[0];
sx q[0];
rz(-1.8183964) q[0];
rz(2.6595751) q[1];
sx q[1];
rz(-2.2299485) q[1];
sx q[1];
rz(2.1673896) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0486748) q[0];
sx q[0];
rz(-2.4705624) q[0];
sx q[0];
rz(-1.2811518) q[0];
rz(2.2902238) q[2];
sx q[2];
rz(-0.15002827) q[2];
sx q[2];
rz(-2.0025557) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43073359) q[1];
sx q[1];
rz(-1.3361592) q[1];
sx q[1];
rz(-2.5905142) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5663381) q[3];
sx q[3];
rz(-1.3674595) q[3];
sx q[3];
rz(-1.4161863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1366068) q[2];
sx q[2];
rz(-2.5544781) q[2];
sx q[2];
rz(2.3919487) q[2];
rz(2.7096115) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(1.3031134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5169446) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(2.5352449) q[0];
rz(-1.8079405) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(2.8820754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57588803) q[0];
sx q[0];
rz(-1.4153) q[0];
sx q[0];
rz(-1.8207654) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9664842) q[2];
sx q[2];
rz(-1.7705743) q[2];
sx q[2];
rz(-1.7410884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9898228) q[1];
sx q[1];
rz(-0.58126106) q[1];
sx q[1];
rz(-0.84110028) q[1];
rz(-2.7068038) q[3];
sx q[3];
rz(-2.090095) q[3];
sx q[3];
rz(-3.0360707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2567265) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(1.5941031) q[2];
rz(-2.3571842) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467658) q[0];
sx q[0];
rz(-1.0972728) q[0];
sx q[0];
rz(2.8821017) q[0];
rz(-1.9208113) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(-1.6993914) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9913015) q[0];
sx q[0];
rz(-2.0388921) q[0];
sx q[0];
rz(0.97458124) q[0];
rz(-pi) q[1];
rz(-0.97564189) q[2];
sx q[2];
rz(-1.4624498) q[2];
sx q[2];
rz(0.25749046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4220548) q[1];
sx q[1];
rz(-1.9214848) q[1];
sx q[1];
rz(-1.3077875) q[1];
rz(-pi) q[2];
rz(1.4846701) q[3];
sx q[3];
rz(-2.4575876) q[3];
sx q[3];
rz(-0.67953426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0356902) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(2.7316459) q[2];
rz(-1.2116872) q[3];
sx q[3];
rz(-0.68710059) q[3];
sx q[3];
rz(1.80779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164417) q[0];
sx q[0];
rz(-0.71617675) q[0];
sx q[0];
rz(0.27405611) q[0];
rz(0.43824276) q[1];
sx q[1];
rz(-1.848105) q[1];
sx q[1];
rz(0.73572198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5925795) q[0];
sx q[0];
rz(-1.4837449) q[0];
sx q[0];
rz(0.75496952) q[0];
x q[1];
rz(-2.3724298) q[2];
sx q[2];
rz(-2.935098) q[2];
sx q[2];
rz(1.002305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3408518) q[1];
sx q[1];
rz(-1.2781118) q[1];
sx q[1];
rz(0.26009788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92145958) q[3];
sx q[3];
rz(-1.9917352) q[3];
sx q[3];
rz(-0.73900797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.938544) q[2];
sx q[2];
rz(-1.6837348) q[2];
sx q[2];
rz(0.61384821) q[2];
rz(-2.7326873) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(2.3992505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78855377) q[0];
sx q[0];
rz(-0.63816324) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(-1.4452112) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(0.17527418) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1244558) q[0];
sx q[0];
rz(-1.5958028) q[0];
sx q[0];
rz(1.7437922) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66368033) q[2];
sx q[2];
rz(-2.366407) q[2];
sx q[2];
rz(-0.41942393) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3172463) q[1];
sx q[1];
rz(-0.69686705) q[1];
sx q[1];
rz(-0.70965712) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95703362) q[3];
sx q[3];
rz(-1.0816469) q[3];
sx q[3];
rz(0.027994284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.96482977) q[2];
sx q[2];
rz(-1.333933) q[2];
sx q[2];
rz(2.1691587) q[2];
rz(-1.6644299) q[3];
sx q[3];
rz(-1.1494613) q[3];
sx q[3];
rz(2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85457388) q[0];
sx q[0];
rz(-3.119097) q[0];
sx q[0];
rz(2.7884685) q[0];
rz(2.1137721) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(-1.0110528) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42732692) q[0];
sx q[0];
rz(-1.9280701) q[0];
sx q[0];
rz(-1.5861804) q[0];
rz(-2.249116) q[2];
sx q[2];
rz(-0.8304285) q[2];
sx q[2];
rz(0.76039808) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87349975) q[1];
sx q[1];
rz(-0.38485369) q[1];
sx q[1];
rz(1.1420239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8854463) q[3];
sx q[3];
rz(-1.2928179) q[3];
sx q[3];
rz(1.6772456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8280243) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(-1.774452) q[2];
rz(1.2683055) q[3];
sx q[3];
rz(-1.4788078) q[3];
sx q[3];
rz(-2.2936599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0276412) q[0];
sx q[0];
rz(-0.60751644) q[0];
sx q[0];
rz(-1.129958) q[0];
rz(2.9226411) q[1];
sx q[1];
rz(-1.4066701) q[1];
sx q[1];
rz(0.91046441) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.509678) q[0];
sx q[0];
rz(-2.4408337) q[0];
sx q[0];
rz(2.6118082) q[0];
rz(-1.4850281) q[2];
sx q[2];
rz(-1.0940486) q[2];
sx q[2];
rz(1.3439182) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95751244) q[1];
sx q[1];
rz(-2.0801447) q[1];
sx q[1];
rz(0.30708509) q[1];
rz(-1.1651785) q[3];
sx q[3];
rz(-2.2542103) q[3];
sx q[3];
rz(0.77714099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3160481) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(-0.82553378) q[2];
rz(2.6801706) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(2.6958444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.5017186) q[0];
sx q[0];
rz(-1.3383144) q[0];
sx q[0];
rz(0.83183944) q[0];
rz(0.72744751) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(0.096253455) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.03503) q[0];
sx q[0];
rz(-2.0756531) q[0];
sx q[0];
rz(-2.6120606) q[0];
rz(-pi) q[1];
rz(2.8873994) q[2];
sx q[2];
rz(-1.6302178) q[2];
sx q[2];
rz(-0.78320247) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2029496) q[1];
sx q[1];
rz(-1.4803107) q[1];
sx q[1];
rz(-1.3590779) q[1];
x q[2];
rz(2.9828137) q[3];
sx q[3];
rz(-1.3769866) q[3];
sx q[3];
rz(-2.4474582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7353797) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(-1.5489138) q[2];
rz(2.5088572) q[3];
sx q[3];
rz(-1.9097208) q[3];
sx q[3];
rz(2.8922141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7525472) q[0];
sx q[0];
rz(-2.1010375) q[0];
sx q[0];
rz(-2.7885875) q[0];
rz(-2.8189335) q[1];
sx q[1];
rz(-0.32538515) q[1];
sx q[1];
rz(-0.37193146) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8648099) q[0];
sx q[0];
rz(-1.7686592) q[0];
sx q[0];
rz(-1.8772621) q[0];
rz(2.6549321) q[2];
sx q[2];
rz(-0.36937215) q[2];
sx q[2];
rz(-2.3104359) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5185753) q[1];
sx q[1];
rz(-1.2378344) q[1];
sx q[1];
rz(-1.740231) q[1];
rz(-1.6149893) q[3];
sx q[3];
rz(-0.72504504) q[3];
sx q[3];
rz(2.8738662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.17452621) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(-1.8444427) q[2];
rz(2.2938812) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(2.1367836) q[3];
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
rz(0.90947718) q[0];
sx q[0];
rz(-1.3790601) q[0];
sx q[0];
rz(-2.7098304) q[0];
rz(-1.3879981) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(0.8193916) q[2];
sx q[2];
rz(-2.1234305) q[2];
sx q[2];
rz(2.7101868) q[2];
rz(1.150849) q[3];
sx q[3];
rz(-2.3227878) q[3];
sx q[3];
rz(-0.25521758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
