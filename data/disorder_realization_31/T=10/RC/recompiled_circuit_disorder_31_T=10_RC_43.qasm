OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52656093) q[0];
sx q[0];
rz(-2.5685413) q[0];
sx q[0];
rz(-0.84258643) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(-1.6834747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18544491) q[0];
sx q[0];
rz(-1.3584134) q[0];
sx q[0];
rz(2.3953715) q[0];
rz(-pi) q[1];
rz(1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(-0.27054271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9021437) q[1];
sx q[1];
rz(-1.4414756) q[1];
sx q[1];
rz(-0.37954482) q[1];
x q[2];
rz(1.6928715) q[3];
sx q[3];
rz(-2.1616518) q[3];
sx q[3];
rz(-1.6216244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52790102) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801385) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.9869841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017529537) q[0];
sx q[0];
rz(-1.5520099) q[0];
sx q[0];
rz(1.5879052) q[0];
x q[1];
rz(-2.7484659) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(1.7413505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0305811) q[1];
sx q[1];
rz(-1.6554553) q[1];
sx q[1];
rz(2.3762523) q[1];
rz(-1.8335908) q[3];
sx q[3];
rz(-1.7495219) q[3];
sx q[3];
rz(0.59584111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(-2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(2.0498958) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0059144817) q[0];
sx q[0];
rz(-1.0251097) q[0];
sx q[0];
rz(-0.90555993) q[0];
rz(-0.91471471) q[2];
sx q[2];
rz(-1.2351742) q[2];
sx q[2];
rz(-0.4682954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6669238) q[1];
sx q[1];
rz(-2.3033934) q[1];
sx q[1];
rz(-1.0679507) q[1];
x q[2];
rz(2.6767119) q[3];
sx q[3];
rz(-2.998623) q[3];
sx q[3];
rz(0.79899341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(-0.88095218) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(2.8312347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096075637) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(1.4287352) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68508673) q[2];
sx q[2];
rz(-1.473046) q[2];
sx q[2];
rz(-0.098066559) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6314108) q[1];
sx q[1];
rz(-1.5893755) q[1];
sx q[1];
rz(-1.9287964) q[1];
rz(-pi) q[2];
rz(-1.6944828) q[3];
sx q[3];
rz(-1.5354772) q[3];
sx q[3];
rz(2.8992821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(3.0854026) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(0.24965832) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10660431) q[0];
sx q[0];
rz(-0.37811324) q[0];
sx q[0];
rz(0.58017054) q[0];
rz(2.9179847) q[2];
sx q[2];
rz(-2.4325271) q[2];
sx q[2];
rz(-2.2774334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.732547) q[1];
sx q[1];
rz(-0.88102075) q[1];
sx q[1];
rz(-1.9802666) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0088483) q[3];
sx q[3];
rz(-1.4268488) q[3];
sx q[3];
rz(2.2465003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(2.8379748) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(2.8836024) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0598037) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(-2.9102737) q[0];
rz(-0.67955534) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(2.213775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7266453) q[1];
sx q[1];
rz(-0.84268314) q[1];
sx q[1];
rz(-0.77264087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8317354) q[3];
sx q[3];
rz(-2.76537) q[3];
sx q[3];
rz(0.46686831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(3.0498665) q[2];
rz(0.84364676) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(2.2475524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-0.86589083) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(3.1076028) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-2.3410428) q[0];
sx q[0];
rz(0.22649015) q[0];
rz(-pi) q[1];
rz(-2.6007973) q[2];
sx q[2];
rz(-2.135709) q[2];
sx q[2];
rz(-2.5164547) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9719203) q[1];
sx q[1];
rz(-1.5202513) q[1];
sx q[1];
rz(2.2488942) q[1];
rz(-2.9339318) q[3];
sx q[3];
rz(-2.5163979) q[3];
sx q[3];
rz(-3.0561662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.6961018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5515585) q[0];
sx q[0];
rz(-0.43338767) q[0];
sx q[0];
rz(1.5584459) q[0];
rz(-pi) q[1];
rz(0.58480279) q[2];
sx q[2];
rz(-1.0324761) q[2];
sx q[2];
rz(-2.1540097) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1576924) q[1];
sx q[1];
rz(-1.0219814) q[1];
sx q[1];
rz(-1.7556612) q[1];
rz(0.94533841) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(-0.64627796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(-1.5173222) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.2517713) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(0.30977419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1447434) q[0];
sx q[0];
rz(-2.092917) q[0];
sx q[0];
rz(-1.7168619) q[0];
rz(1.6236213) q[2];
sx q[2];
rz(-0.17715684) q[2];
sx q[2];
rz(-1.0747386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5619547) q[1];
sx q[1];
rz(-1.5127752) q[1];
sx q[1];
rz(-1.3647563) q[1];
x q[2];
rz(-1.8272607) q[3];
sx q[3];
rz(-2.1520352) q[3];
sx q[3];
rz(-2.6687711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.4996128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033894) q[0];
sx q[0];
rz(-1.8614385) q[0];
sx q[0];
rz(-3.035726) q[0];
rz(-pi) q[1];
rz(-1.502938) q[2];
sx q[2];
rz(-0.80791622) q[2];
sx q[2];
rz(2.9615336) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2227576) q[1];
sx q[1];
rz(-1.3593874) q[1];
sx q[1];
rz(-1.097015) q[1];
x q[2];
rz(-1.3979982) q[3];
sx q[3];
rz(-1.5683335) q[3];
sx q[3];
rz(2.5415004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(2.3161841) q[2];
sx q[2];
rz(-2.5054629) q[2];
sx q[2];
rz(2.7873743) q[2];
rz(2.0675038) q[3];
sx q[3];
rz(-0.91377331) q[3];
sx q[3];
rz(2.5627315) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
