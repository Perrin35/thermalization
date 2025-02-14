OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(-1.2824143) q[0];
sx q[0];
rz(-0.13089827) q[0];
rz(-2.6137597) q[1];
sx q[1];
rz(-0.37017828) q[1];
sx q[1];
rz(0.17732492) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89761868) q[0];
sx q[0];
rz(-1.8332687) q[0];
sx q[0];
rz(-0.7926244) q[0];
rz(-1.2700291) q[2];
sx q[2];
rz(-1.352515) q[2];
sx q[2];
rz(-2.0576365) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51715358) q[1];
sx q[1];
rz(-0.53688795) q[1];
sx q[1];
rz(-1.0978903) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8683441) q[3];
sx q[3];
rz(-2.1299703) q[3];
sx q[3];
rz(2.3756053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6191972) q[2];
sx q[2];
rz(-2.0530687) q[2];
sx q[2];
rz(1.8001451) q[2];
rz(1.7893192) q[3];
sx q[3];
rz(-1.5315346) q[3];
sx q[3];
rz(1.8488098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.751048) q[0];
sx q[0];
rz(-1.229267) q[0];
sx q[0];
rz(2.6265662) q[0];
rz(2.7797049) q[1];
sx q[1];
rz(-1.7444976) q[1];
sx q[1];
rz(-0.99641689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85595713) q[0];
sx q[0];
rz(-2.8295316) q[0];
sx q[0];
rz(0.86459222) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9662232) q[2];
sx q[2];
rz(-1.9992454) q[2];
sx q[2];
rz(-2.1074949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39033088) q[1];
sx q[1];
rz(-1.0511025) q[1];
sx q[1];
rz(3.0355176) q[1];
x q[2];
rz(2.6985136) q[3];
sx q[3];
rz(-2.0650117) q[3];
sx q[3];
rz(-2.6576603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4126052) q[2];
sx q[2];
rz(-0.94634405) q[2];
sx q[2];
rz(1.4808572) q[2];
rz(2.6442773) q[3];
sx q[3];
rz(-2.1931084) q[3];
sx q[3];
rz(-0.724154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5363252) q[0];
sx q[0];
rz(-1.241642) q[0];
sx q[0];
rz(1.190825) q[0];
rz(-0.39769998) q[1];
sx q[1];
rz(-1.560248) q[1];
sx q[1];
rz(-2.2427028) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1440292) q[0];
sx q[0];
rz(-1.231651) q[0];
sx q[0];
rz(-3.0710286) q[0];
rz(-2.7238893) q[2];
sx q[2];
rz(-2.933393) q[2];
sx q[2];
rz(2.4793159) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6582322) q[1];
sx q[1];
rz(-1.0387324) q[1];
sx q[1];
rz(3.0897555) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0895899) q[3];
sx q[3];
rz(-1.7214603) q[3];
sx q[3];
rz(2.0637531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1206104) q[2];
sx q[2];
rz(-1.9526498) q[2];
sx q[2];
rz(2.9494542) q[2];
rz(-0.7297248) q[3];
sx q[3];
rz(-0.14484043) q[3];
sx q[3];
rz(-2.5016968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14908734) q[0];
sx q[0];
rz(-0.75228107) q[0];
sx q[0];
rz(-1.8335023) q[0];
rz(-2.2970842) q[1];
sx q[1];
rz(-1.4627855) q[1];
sx q[1];
rz(0.62215296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7956411) q[0];
sx q[0];
rz(-1.9824195) q[0];
sx q[0];
rz(2.6108885) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18537955) q[2];
sx q[2];
rz(-1.6246008) q[2];
sx q[2];
rz(2.0410694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9355008) q[1];
sx q[1];
rz(-2.969944) q[1];
sx q[1];
rz(2.8354581) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34028168) q[3];
sx q[3];
rz(-0.33512402) q[3];
sx q[3];
rz(1.7652896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1239803) q[2];
sx q[2];
rz(-1.3717185) q[2];
sx q[2];
rz(-2.2464216) q[2];
rz(-0.34919843) q[3];
sx q[3];
rz(-2.5694191) q[3];
sx q[3];
rz(-2.9077742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2772086) q[0];
sx q[0];
rz(-1.9338436) q[0];
sx q[0];
rz(2.5982017) q[0];
rz(3.0282989) q[1];
sx q[1];
rz(-1.7430867) q[1];
sx q[1];
rz(-1.2921804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3709576) q[0];
sx q[0];
rz(-1.1877302) q[0];
sx q[0];
rz(-2.2983453) q[0];
rz(-0.92083709) q[2];
sx q[2];
rz(-1.2652768) q[2];
sx q[2];
rz(0.68543514) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1317392) q[1];
sx q[1];
rz(-0.92127548) q[1];
sx q[1];
rz(-2.3765537) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4025147) q[3];
sx q[3];
rz(-1.0024286) q[3];
sx q[3];
rz(3.0773602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8285404) q[2];
sx q[2];
rz(-1.3581759) q[2];
sx q[2];
rz(-2.6521315) q[2];
rz(2.475907) q[3];
sx q[3];
rz(-0.61167115) q[3];
sx q[3];
rz(-0.78589511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2419484) q[0];
sx q[0];
rz(-0.90270942) q[0];
sx q[0];
rz(-3.0808501) q[0];
rz(1.863106) q[1];
sx q[1];
rz(-2.2272031) q[1];
sx q[1];
rz(-1.0260822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4961012) q[0];
sx q[0];
rz(-2.0466986) q[0];
sx q[0];
rz(0.93604139) q[0];
rz(-1.0025239) q[2];
sx q[2];
rz(-1.2934877) q[2];
sx q[2];
rz(-0.31686767) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85426165) q[1];
sx q[1];
rz(-1.8214487) q[1];
sx q[1];
rz(-1.7118368) q[1];
rz(-2.6164651) q[3];
sx q[3];
rz(-2.5487102) q[3];
sx q[3];
rz(-1.2660932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25980276) q[2];
sx q[2];
rz(-1.5851574) q[2];
sx q[2];
rz(-3.0628824) q[2];
rz(-1.4824661) q[3];
sx q[3];
rz(-0.31646287) q[3];
sx q[3];
rz(-0.44221529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92597961) q[0];
sx q[0];
rz(-0.38023606) q[0];
sx q[0];
rz(-1.6967787) q[0];
rz(0.80698693) q[1];
sx q[1];
rz(-1.746256) q[1];
sx q[1];
rz(3.1296465) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1211469) q[0];
sx q[0];
rz(-0.63948005) q[0];
sx q[0];
rz(1.5259652) q[0];
rz(-pi) q[1];
rz(-0.50759727) q[2];
sx q[2];
rz(-1.6656808) q[2];
sx q[2];
rz(-2.0213057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.641549) q[1];
sx q[1];
rz(-2.1807359) q[1];
sx q[1];
rz(-0.088492101) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1720522) q[3];
sx q[3];
rz(-1.4280025) q[3];
sx q[3];
rz(-0.90313426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.547946) q[2];
sx q[2];
rz(-2.839489) q[2];
sx q[2];
rz(2.3054403) q[2];
rz(-3.0892843) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(-0.93594319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.1702105) q[0];
sx q[0];
rz(-2.2844071) q[0];
sx q[0];
rz(0.49825391) q[0];
rz(1.0055297) q[1];
sx q[1];
rz(-1.6906831) q[1];
sx q[1];
rz(3.0301869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69499066) q[0];
sx q[0];
rz(-1.989762) q[0];
sx q[0];
rz(-0.10563157) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6314503) q[2];
sx q[2];
rz(-2.6594793) q[2];
sx q[2];
rz(3.1010266) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3058779) q[1];
sx q[1];
rz(-1.8566161) q[1];
sx q[1];
rz(2.3217683) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8008055) q[3];
sx q[3];
rz(-0.36431815) q[3];
sx q[3];
rz(0.43393578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19693836) q[2];
sx q[2];
rz(-1.3613746) q[2];
sx q[2];
rz(1.4062175) q[2];
rz(-1.9841291) q[3];
sx q[3];
rz(-0.99298733) q[3];
sx q[3];
rz(-1.5520613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.120753) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(1.7027759) q[0];
rz(2.2976047) q[1];
sx q[1];
rz(-1.3071209) q[1];
sx q[1];
rz(0.2058952) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5000847) q[0];
sx q[0];
rz(-2.5091268) q[0];
sx q[0];
rz(-3.0728805) q[0];
x q[1];
rz(1.6240142) q[2];
sx q[2];
rz(-0.88376617) q[2];
sx q[2];
rz(1.991697) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2515427) q[1];
sx q[1];
rz(-1.6487507) q[1];
sx q[1];
rz(-0.13549094) q[1];
rz(-2.467942) q[3];
sx q[3];
rz(-0.60551548) q[3];
sx q[3];
rz(-0.57693627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0042808) q[2];
sx q[2];
rz(-2.0840804) q[2];
sx q[2];
rz(1.6005969) q[2];
rz(-0.48504034) q[3];
sx q[3];
rz(-1.78616) q[3];
sx q[3];
rz(-3.0276827) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4892905) q[0];
sx q[0];
rz(-0.052209608) q[0];
sx q[0];
rz(1.2963649) q[0];
rz(-1.0507874) q[1];
sx q[1];
rz(-1.4605582) q[1];
sx q[1];
rz(0.60660648) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8826128) q[0];
sx q[0];
rz(-1.4073696) q[0];
sx q[0];
rz(0.033074065) q[0];
x q[1];
rz(-1.6982618) q[2];
sx q[2];
rz(-2.1268788) q[2];
sx q[2];
rz(-3.0024204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21571017) q[1];
sx q[1];
rz(-1.6318738) q[1];
sx q[1];
rz(-0.30558773) q[1];
rz(-0.38207558) q[3];
sx q[3];
rz(-0.30508546) q[3];
sx q[3];
rz(0.24972734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4198833) q[2];
sx q[2];
rz(-1.1315283) q[2];
sx q[2];
rz(0.69880542) q[2];
rz(-1.6803668) q[3];
sx q[3];
rz(-1.0137839) q[3];
sx q[3];
rz(-0.47718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2496495) q[0];
sx q[0];
rz(-1.5105381) q[0];
sx q[0];
rz(1.426209) q[0];
rz(-2.7632948) q[1];
sx q[1];
rz(-2.5143647) q[1];
sx q[1];
rz(-0.73077269) q[1];
rz(-1.0155914) q[2];
sx q[2];
rz(-2.5980661) q[2];
sx q[2];
rz(-0.51474434) q[2];
rz(1.6017687) q[3];
sx q[3];
rz(-2.1270558) q[3];
sx q[3];
rz(1.1598641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
