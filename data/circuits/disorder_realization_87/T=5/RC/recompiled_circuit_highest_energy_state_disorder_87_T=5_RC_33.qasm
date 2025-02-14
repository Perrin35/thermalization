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
rz(5.7558007) q[0];
sx q[0];
rz(5.7100073) q[0];
sx q[0];
rz(10.041458) q[0];
rz(2.1332027) q[1];
sx q[1];
rz(-0.73928666) q[1];
sx q[1];
rz(-2.8032805) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9341921) q[0];
sx q[0];
rz(-1.8131885) q[0];
sx q[0];
rz(-1.6079748) q[0];
rz(-2.7213547) q[2];
sx q[2];
rz(-2.3048688) q[2];
sx q[2];
rz(-0.24009304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9788614) q[1];
sx q[1];
rz(-0.56038364) q[1];
sx q[1];
rz(-0.25516971) q[1];
x q[2];
rz(-0.81387384) q[3];
sx q[3];
rz(-1.4424837) q[3];
sx q[3];
rz(-1.8112884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7484635) q[2];
sx q[2];
rz(-1.7270361) q[2];
sx q[2];
rz(0.65790042) q[2];
rz(-2.6864478) q[3];
sx q[3];
rz(-2.9908266) q[3];
sx q[3];
rz(1.3078825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748134) q[0];
sx q[0];
rz(-1.4115189) q[0];
sx q[0];
rz(0.28208062) q[0];
rz(-0.2433978) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(-0.74877053) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.153745) q[0];
sx q[0];
rz(-1.3778338) q[0];
sx q[0];
rz(0.70867507) q[0];
x q[1];
rz(1.6924627) q[2];
sx q[2];
rz(-1.5632544) q[2];
sx q[2];
rz(-2.6040524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2255324) q[1];
sx q[1];
rz(-2.6301363) q[1];
sx q[1];
rz(1.3813853) q[1];
rz(-pi) q[2];
rz(0.33269675) q[3];
sx q[3];
rz(-2.930958) q[3];
sx q[3];
rz(-3.0300107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8613209) q[2];
sx q[2];
rz(-1.5254285) q[2];
sx q[2];
rz(-1.8005499) q[2];
rz(1.9848112) q[3];
sx q[3];
rz(-0.78514922) q[3];
sx q[3];
rz(-2.3780499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083104) q[0];
sx q[0];
rz(-0.16591993) q[0];
sx q[0];
rz(-2.4842343) q[0];
rz(-1.1454469) q[1];
sx q[1];
rz(-0.95445389) q[1];
sx q[1];
rz(-2.3232536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2812202) q[0];
sx q[0];
rz(-1.91433) q[0];
sx q[0];
rz(-1.3168797) q[0];
rz(1.255515) q[2];
sx q[2];
rz(-1.6161801) q[2];
sx q[2];
rz(-0.78901828) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1166621) q[1];
sx q[1];
rz(-2.0567472) q[1];
sx q[1];
rz(2.3930753) q[1];
x q[2];
rz(3.0679296) q[3];
sx q[3];
rz(-2.0233002) q[3];
sx q[3];
rz(1.0301925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1546617) q[2];
sx q[2];
rz(-1.8795452) q[2];
sx q[2];
rz(-1.3531125) q[2];
rz(-2.1766369) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(-2.7195948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6956536) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(-1.2009784) q[0];
rz(-3.1383842) q[1];
sx q[1];
rz(-1.4326347) q[1];
sx q[1];
rz(-2.9046955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85962109) q[0];
sx q[0];
rz(-1.0399204) q[0];
sx q[0];
rz(2.6174699) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9658674) q[2];
sx q[2];
rz(-0.64470664) q[2];
sx q[2];
rz(0.75632655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5665) q[1];
sx q[1];
rz(-1.3761576) q[1];
sx q[1];
rz(-0.95928488) q[1];
x q[2];
rz(-0.77787106) q[3];
sx q[3];
rz(-0.98613769) q[3];
sx q[3];
rz(-0.28077048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0420456) q[2];
sx q[2];
rz(-1.2914265) q[2];
sx q[2];
rz(-2.715204) q[2];
rz(-2.1060627) q[3];
sx q[3];
rz(-1.6798881) q[3];
sx q[3];
rz(2.5199913) q[3];
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
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7972888) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(0.34991831) q[0];
rz(1.2087076) q[1];
sx q[1];
rz(-1.2827001) q[1];
sx q[1];
rz(-3.1399609) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37536538) q[0];
sx q[0];
rz(-1.3266801) q[0];
sx q[0];
rz(-1.7422416) q[0];
rz(-0.12589215) q[2];
sx q[2];
rz(-2.841904) q[2];
sx q[2];
rz(-2.3131049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.092068191) q[1];
sx q[1];
rz(-2.9584269) q[1];
sx q[1];
rz(1.8335672) q[1];
x q[2];
rz(-2.9309421) q[3];
sx q[3];
rz(-1.8351166) q[3];
sx q[3];
rz(1.5624865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1899167) q[2];
sx q[2];
rz(-2.2168171) q[2];
sx q[2];
rz(0.16119371) q[2];
rz(-0.4392043) q[3];
sx q[3];
rz(-2.3930211) q[3];
sx q[3];
rz(-2.2883033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.31140232) q[0];
sx q[0];
rz(-1.7364194) q[0];
sx q[0];
rz(2.3468974) q[0];
rz(-1.3658124) q[1];
sx q[1];
rz(-0.8005442) q[1];
sx q[1];
rz(-2.6672003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77664033) q[0];
sx q[0];
rz(-2.4863003) q[0];
sx q[0];
rz(2.220015) q[0];
x q[1];
rz(2.5854255) q[2];
sx q[2];
rz(-0.63236299) q[2];
sx q[2];
rz(1.7025422) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77955708) q[1];
sx q[1];
rz(-1.0048702) q[1];
sx q[1];
rz(1.0841682) q[1];
rz(-pi) q[2];
rz(-0.9871363) q[3];
sx q[3];
rz(-1.8963834) q[3];
sx q[3];
rz(-1.1914355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5279493) q[2];
sx q[2];
rz(-1.8972634) q[2];
sx q[2];
rz(2.6141686) q[2];
rz(-2.534965) q[3];
sx q[3];
rz(-2.9546723) q[3];
sx q[3];
rz(-1.0724148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793176) q[0];
sx q[0];
rz(-2.9476705) q[0];
sx q[0];
rz(-2.344017) q[0];
rz(-0.89353117) q[1];
sx q[1];
rz(-0.83502665) q[1];
sx q[1];
rz(-0.28018793) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3994795) q[0];
sx q[0];
rz(-0.43599883) q[0];
sx q[0];
rz(-1.9943555) q[0];
rz(-pi) q[1];
rz(-2.3827219) q[2];
sx q[2];
rz(-2.3312097) q[2];
sx q[2];
rz(-2.9019865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60498991) q[1];
sx q[1];
rz(-2.0532673) q[1];
sx q[1];
rz(0.014181344) q[1];
x q[2];
rz(2.7817621) q[3];
sx q[3];
rz(-2.1323279) q[3];
sx q[3];
rz(-2.7240745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4013275) q[2];
sx q[2];
rz(-1.3046616) q[2];
sx q[2];
rz(1.2124445) q[2];
rz(2.9017743) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(2.5734606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577393) q[0];
sx q[0];
rz(-1.4171866) q[0];
sx q[0];
rz(1.8200112) q[0];
rz(-2.9603738) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(-0.063057335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6609333) q[0];
sx q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(0.78137915) q[0];
rz(-1.2853773) q[2];
sx q[2];
rz(-1.4865424) q[2];
sx q[2];
rz(0.38049305) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93062295) q[1];
sx q[1];
rz(-1.2377059) q[1];
sx q[1];
rz(1.3107783) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1386915) q[3];
sx q[3];
rz(-1.408268) q[3];
sx q[3];
rz(1.8123466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41398373) q[2];
sx q[2];
rz(-2.0435645) q[2];
sx q[2];
rz(-1.6723527) q[2];
rz(-1.3937048) q[3];
sx q[3];
rz(-1.4398451) q[3];
sx q[3];
rz(0.20974717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53836981) q[0];
sx q[0];
rz(-0.61036888) q[0];
sx q[0];
rz(2.1446153) q[0];
rz(-0.793055) q[1];
sx q[1];
rz(-1.7231562) q[1];
sx q[1];
rz(-2.0319895) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9859409) q[0];
sx q[0];
rz(-1.8619259) q[0];
sx q[0];
rz(-0.58737866) q[0];
rz(0.72203861) q[2];
sx q[2];
rz(-1.8991611) q[2];
sx q[2];
rz(1.2102845) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3396757) q[1];
sx q[1];
rz(-2.0387421) q[1];
sx q[1];
rz(-0.94593327) q[1];
rz(-pi) q[2];
rz(-1.4023191) q[3];
sx q[3];
rz(-2.2003897) q[3];
sx q[3];
rz(-1.633267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41254607) q[2];
sx q[2];
rz(-1.6763326) q[2];
sx q[2];
rz(-1.868978) q[2];
rz(-1.1254958) q[3];
sx q[3];
rz(-0.19386217) q[3];
sx q[3];
rz(-3.061749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37128714) q[0];
sx q[0];
rz(-1.1530387) q[0];
sx q[0];
rz(-3.0433997) q[0];
rz(1.498361) q[1];
sx q[1];
rz(-1.1474835) q[1];
sx q[1];
rz(0.8383382) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6107008) q[0];
sx q[0];
rz(-0.77741388) q[0];
sx q[0];
rz(-0.25334187) q[0];
rz(0.25213253) q[2];
sx q[2];
rz(-1.5169391) q[2];
sx q[2];
rz(0.60562741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.34217087) q[1];
sx q[1];
rz(-1.2890745) q[1];
sx q[1];
rz(-2.0215624) q[1];
rz(1.3499979) q[3];
sx q[3];
rz(-1.1482557) q[3];
sx q[3];
rz(-2.411946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7958293) q[2];
sx q[2];
rz(-0.886262) q[2];
sx q[2];
rz(0.32540992) q[2];
rz(0.77643967) q[3];
sx q[3];
rz(-2.1263945) q[3];
sx q[3];
rz(-2.5344892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0008739) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(2.8554032) q[1];
sx q[1];
rz(-1.9356526) q[1];
sx q[1];
rz(-2.0331358) q[1];
rz(2.3807965) q[2];
sx q[2];
rz(-1.2789055) q[2];
sx q[2];
rz(0.61550752) q[2];
rz(2.6558032) q[3];
sx q[3];
rz(-2.494907) q[3];
sx q[3];
rz(-3.1068556) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
