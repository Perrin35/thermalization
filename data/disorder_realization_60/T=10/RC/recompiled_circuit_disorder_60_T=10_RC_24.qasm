OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.94937593) q[0];
sx q[0];
rz(-1.047171) q[0];
sx q[0];
rz(0.068724364) q[0];
rz(1.7460495) q[1];
sx q[1];
rz(4.6739251) q[1];
sx q[1];
rz(8.2164017) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8954593) q[0];
sx q[0];
rz(-1.6172505) q[0];
sx q[0];
rz(-1.6669271) q[0];
rz(-pi) q[1];
rz(3.0460998) q[2];
sx q[2];
rz(-3.0152233) q[2];
sx q[2];
rz(0.40590826) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2627416) q[1];
sx q[1];
rz(-0.95572119) q[1];
sx q[1];
rz(1.4351074) q[1];
rz(-pi) q[2];
rz(0.73322202) q[3];
sx q[3];
rz(-0.65922046) q[3];
sx q[3];
rz(-2.1472907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8157114) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(-0.69774929) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(0.74716032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.117347) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(-1.1741937) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(0.29719621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8288119) q[0];
sx q[0];
rz(-0.096185616) q[0];
sx q[0];
rz(2.0975153) q[0];
rz(-2.6058795) q[2];
sx q[2];
rz(-0.78768724) q[2];
sx q[2];
rz(2.1802769) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9651523) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(-0.39217197) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3791749) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16195665) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(2.5276108) q[2];
rz(0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(1.0959594) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(0.30763787) q[0];
rz(-2.3930507) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(2.3017853) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52536406) q[0];
sx q[0];
rz(-1.1601686) q[0];
sx q[0];
rz(-1.8812268) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4762127) q[2];
sx q[2];
rz(-1.9087831) q[2];
sx q[2];
rz(-2.4199977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0864598) q[1];
sx q[1];
rz(-1.2659402) q[1];
sx q[1];
rz(2.6677368) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3358467) q[3];
sx q[3];
rz(-2.4953105) q[3];
sx q[3];
rz(0.39138734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(1.8236558) q[2];
rz(-1.2157724) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3777305) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(2.6960301) q[0];
rz(-0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(2.1760118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5641898) q[0];
sx q[0];
rz(-1.2947384) q[0];
sx q[0];
rz(0.77703707) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1104286) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(2.503501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.021918745) q[1];
sx q[1];
rz(-1.5389171) q[1];
sx q[1];
rz(-1.8244513) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2772917) q[3];
sx q[3];
rz(-2.9923277) q[3];
sx q[3];
rz(-0.60993689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(0.92932534) q[2];
rz(2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(1.8575645) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(1.0481542) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832344) q[0];
sx q[0];
rz(-2.2336707) q[0];
sx q[0];
rz(-2.5213084) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5571312) q[2];
sx q[2];
rz(-0.9218773) q[2];
sx q[2];
rz(1.6780168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3874665) q[1];
sx q[1];
rz(-1.4782463) q[1];
sx q[1];
rz(0.32317633) q[1];
x q[2];
rz(0.9976451) q[3];
sx q[3];
rz(-0.76612681) q[3];
sx q[3];
rz(-1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(1.0926931) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(-1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95935217) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(-2.545488) q[0];
rz(1.6456564) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.8966819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0743474) q[0];
sx q[0];
rz(-1.0676358) q[0];
sx q[0];
rz(-2.7748681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55576022) q[2];
sx q[2];
rz(-2.5854977) q[2];
sx q[2];
rz(0.13242002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4391172) q[1];
sx q[1];
rz(-2.2687952) q[1];
sx q[1];
rz(0.67081397) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8935029) q[3];
sx q[3];
rz(-1.7956942) q[3];
sx q[3];
rz(3.0450862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46463075) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(-2.8616469) q[0];
rz(1.6784558) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(0.25269145) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0137018) q[0];
sx q[0];
rz(-1.2780006) q[0];
sx q[0];
rz(3.091759) q[0];
rz(-pi) q[1];
rz(1.5451317) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(0.54021013) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45526866) q[1];
sx q[1];
rz(-0.9372006) q[1];
sx q[1];
rz(0.92647657) q[1];
rz(-pi) q[2];
rz(-2.2884376) q[3];
sx q[3];
rz(-2.1907638) q[3];
sx q[3];
rz(0.48228797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32020405) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(-1.5318711) q[2];
rz(-1.1931217) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(-1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367592) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(1.6554792) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(-2.862646) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3725961) q[0];
sx q[0];
rz(-2.0645752) q[0];
sx q[0];
rz(-2.2934224) q[0];
rz(-pi) q[1];
rz(-2.0729162) q[2];
sx q[2];
rz(-1.4197822) q[2];
sx q[2];
rz(2.6627024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3153509) q[1];
sx q[1];
rz(-0.8289753) q[1];
sx q[1];
rz(2.6226603) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7989743) q[3];
sx q[3];
rz(-2.5573213) q[3];
sx q[3];
rz(-2.1669471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.802861) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779697) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(1.8883702) q[0];
rz(-2.5190917) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(1.1463096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.779594) q[0];
sx q[0];
rz(-2.3906374) q[0];
sx q[0];
rz(-0.10303084) q[0];
rz(2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(-2.7149372) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6625001) q[1];
sx q[1];
rz(-0.88516419) q[1];
sx q[1];
rz(-3.0467266) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54543145) q[3];
sx q[3];
rz(-0.71119961) q[3];
sx q[3];
rz(1.4382854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15381947) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(-2.1257607) q[2];
rz(-0.9097957) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(-0.36809665) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(-2.2968538) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(1.3815809) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68122278) q[0];
sx q[0];
rz(-1.5813706) q[0];
sx q[0];
rz(-1.4282385) q[0];
rz(-pi) q[1];
rz(-0.064247473) q[2];
sx q[2];
rz(-0.99284961) q[2];
sx q[2];
rz(0.77677514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(1.0113869) q[1];
x q[2];
rz(1.3515477) q[3];
sx q[3];
rz(-1.0478813) q[3];
sx q[3];
rz(3.074309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0766729) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(-2.1515965) q[2];
rz(-0.89407095) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(-2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0513231) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(1.3399711) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(0.71628911) q[2];
sx q[2];
rz(-1.204797) q[2];
sx q[2];
rz(2.7439678) q[2];
rz(-2.7064825) q[3];
sx q[3];
rz(-2.1463487) q[3];
sx q[3];
rz(-1.6456732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];