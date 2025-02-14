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
rz(1.4229245) q[0];
sx q[0];
rz(-2.0473502) q[0];
sx q[0];
rz(-0.25804582) q[0];
rz(-0.99335042) q[1];
sx q[1];
rz(-1.8408096) q[1];
sx q[1];
rz(2.8859477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48707552) q[0];
sx q[0];
rz(-2.1571113) q[0];
sx q[0];
rz(-2.4357027) q[0];
rz(-0.3886224) q[2];
sx q[2];
rz(-1.6291233) q[2];
sx q[2];
rz(1.3921392) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.930525) q[1];
sx q[1];
rz(-1.3881359) q[1];
sx q[1];
rz(-0.069932368) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3575451) q[3];
sx q[3];
rz(-2.0590326) q[3];
sx q[3];
rz(2.8247716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(1.1154491) q[2];
rz(3.1274146) q[3];
sx q[3];
rz(-1.7970128) q[3];
sx q[3];
rz(0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2071335) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(1.8943262) q[0];
rz(0.44218749) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(-0.8173379) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8431664) q[0];
sx q[0];
rz(-2.0930556) q[0];
sx q[0];
rz(0.22163117) q[0];
rz(-2.1139718) q[2];
sx q[2];
rz(-2.6386542) q[2];
sx q[2];
rz(-2.3523112) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0603588) q[1];
sx q[1];
rz(-2.4801755) q[1];
sx q[1];
rz(0.27137406) q[1];
rz(2.2226187) q[3];
sx q[3];
rz(-0.92286829) q[3];
sx q[3];
rz(-2.1369262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77357972) q[2];
sx q[2];
rz(-0.37675884) q[2];
sx q[2];
rz(-2.9212941) q[2];
rz(-1.9593272) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(3.0784472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050506266) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(0.41172045) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(-2.128111) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6665812) q[0];
sx q[0];
rz(-1.555976) q[0];
sx q[0];
rz(-2.8885452) q[0];
x q[1];
rz(0.13229741) q[2];
sx q[2];
rz(-1.1331285) q[2];
sx q[2];
rz(-2.3486111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.93597465) q[1];
sx q[1];
rz(-2.0009319) q[1];
sx q[1];
rz(0.79367743) q[1];
x q[2];
rz(0.8441505) q[3];
sx q[3];
rz(-1.6682079) q[3];
sx q[3];
rz(1.0507492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1608405) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(-1.8864924) q[2];
rz(0.27215019) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(-3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.960152) q[0];
sx q[0];
rz(-2.20521) q[0];
sx q[0];
rz(-0.61035672) q[0];
rz(-2.4329674) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(1.5707387) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40442586) q[0];
sx q[0];
rz(-2.9350087) q[0];
sx q[0];
rz(-2.4086359) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4518634) q[2];
sx q[2];
rz(-1.2858675) q[2];
sx q[2];
rz(0.97078427) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8766986) q[1];
sx q[1];
rz(-1.2994517) q[1];
sx q[1];
rz(1.62074) q[1];
x q[2];
rz(0.80067486) q[3];
sx q[3];
rz(-1.8599556) q[3];
sx q[3];
rz(1.9137933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9307956) q[2];
sx q[2];
rz(-2.1288629) q[2];
sx q[2];
rz(-0.55541682) q[2];
rz(-2.4173315) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50437462) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(2.7727238) q[0];
rz(-2.8952307) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(-1.4341644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61489366) q[0];
sx q[0];
rz(-1.5656359) q[0];
sx q[0];
rz(2.3847488) q[0];
x q[1];
rz(-1.3952012) q[2];
sx q[2];
rz(-2.0101974) q[2];
sx q[2];
rz(-2.0134157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7692657) q[1];
sx q[1];
rz(-0.40180909) q[1];
sx q[1];
rz(0.010784464) q[1];
rz(-0.09999545) q[3];
sx q[3];
rz(-1.7971562) q[3];
sx q[3];
rz(2.5842427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.65115923) q[2];
sx q[2];
rz(-1.3110524) q[2];
sx q[2];
rz(1.0820214) q[2];
rz(-2.3467482) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(1.883551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.865888) q[0];
sx q[0];
rz(-3.0401433) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(-2.1547735) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(2.2056244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075386062) q[0];
sx q[0];
rz(-2.6200326) q[0];
sx q[0];
rz(-2.2564476) q[0];
x q[1];
rz(2.7924839) q[2];
sx q[2];
rz(-2.287068) q[2];
sx q[2];
rz(0.30308613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0430345) q[1];
sx q[1];
rz(-1.9717311) q[1];
sx q[1];
rz(-0.75717302) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6020847) q[3];
sx q[3];
rz(-2.5415321) q[3];
sx q[3];
rz(0.72573435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.662107) q[2];
sx q[2];
rz(-2.0027436) q[2];
sx q[2];
rz(-0.34995079) q[2];
rz(0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(2.2902655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418942) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(-0.30174524) q[0];
rz(-0.31983495) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(2.005827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1162735) q[0];
sx q[0];
rz(-0.92705446) q[0];
sx q[0];
rz(-1.5035065) q[0];
rz(0.91355027) q[2];
sx q[2];
rz(-0.74591178) q[2];
sx q[2];
rz(-0.28829703) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37138501) q[1];
sx q[1];
rz(-1.6504297) q[1];
sx q[1];
rz(-0.97038986) q[1];
rz(-pi) q[2];
rz(-2.3507422) q[3];
sx q[3];
rz(-1.491427) q[3];
sx q[3];
rz(-0.96109238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4096421) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(0.14275924) q[2];
rz(1.7536633) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(-0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9614354) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(-2.5855682) q[0];
rz(-0.22932886) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(-1.75846) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5189603) q[0];
sx q[0];
rz(-0.83370249) q[0];
sx q[0];
rz(1.5234768) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6298619) q[2];
sx q[2];
rz(-1.1129654) q[2];
sx q[2];
rz(-0.77564592) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38249967) q[1];
sx q[1];
rz(-1.1611995) q[1];
sx q[1];
rz(0.31633693) q[1];
x q[2];
rz(-2.4501958) q[3];
sx q[3];
rz(-2.5775839) q[3];
sx q[3];
rz(-2.9897652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6250299) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(1.2410835) q[2];
rz(-0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(-0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1271707) q[0];
sx q[0];
rz(-1.4143455) q[0];
sx q[0];
rz(0.14778368) q[0];
rz(0.5087018) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(-2.7630189) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157279) q[0];
sx q[0];
rz(-1.8047389) q[0];
sx q[0];
rz(-0.42629231) q[0];
rz(0.97855391) q[2];
sx q[2];
rz(-1.2662953) q[2];
sx q[2];
rz(-1.026012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3264824) q[1];
sx q[1];
rz(-1.8811783) q[1];
sx q[1];
rz(1.2182477) q[1];
rz(-pi) q[2];
x q[2];
rz(2.380563) q[3];
sx q[3];
rz(-1.4806403) q[3];
sx q[3];
rz(0.48616274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1754237) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(1.8625205) q[2];
rz(1.7133948) q[3];
sx q[3];
rz(-1.0019852) q[3];
sx q[3];
rz(0.14555791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(1.3928423) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(0.064099126) q[0];
rz(-2.6129258) q[1];
sx q[1];
rz(-1.8730947) q[1];
sx q[1];
rz(2.166523) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9720358) q[0];
sx q[0];
rz(-1.0825048) q[0];
sx q[0];
rz(1.2704865) q[0];
rz(2.5717402) q[2];
sx q[2];
rz(-2.0790711) q[2];
sx q[2];
rz(-3.0076671) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1499612) q[1];
sx q[1];
rz(-1.35222) q[1];
sx q[1];
rz(-2.6944955) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7998135) q[3];
sx q[3];
rz(-1.8088565) q[3];
sx q[3];
rz(2.2464667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9958682) q[2];
sx q[2];
rz(-2.4805562) q[2];
sx q[2];
rz(-2.5346942) q[2];
rz(2.8912344) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(-0.50601602) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.09457) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(2.6429214) q[1];
sx q[1];
rz(-1.5892727) q[1];
sx q[1];
rz(-1.7472063) q[1];
rz(-2.5468536) q[2];
sx q[2];
rz(-0.95437106) q[2];
sx q[2];
rz(0.30827733) q[2];
rz(-1.7767033) q[3];
sx q[3];
rz(-1.4510703) q[3];
sx q[3];
rz(2.5616796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
