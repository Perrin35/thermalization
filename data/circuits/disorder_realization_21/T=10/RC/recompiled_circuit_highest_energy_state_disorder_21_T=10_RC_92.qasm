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
rz(0.30581185) q[0];
sx q[0];
rz(-2.5078791) q[0];
sx q[0];
rz(0.60854882) q[0];
rz(2.4701056) q[1];
sx q[1];
rz(-0.66317135) q[1];
sx q[1];
rz(1.6821678) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6240294) q[0];
sx q[0];
rz(-1.3362954) q[0];
sx q[0];
rz(1.0453216) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9179108) q[2];
sx q[2];
rz(-2.0111901) q[2];
sx q[2];
rz(0.5612095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0406831) q[1];
sx q[1];
rz(-0.95129943) q[1];
sx q[1];
rz(-2.6454974) q[1];
x q[2];
rz(2.9608034) q[3];
sx q[3];
rz(-0.76813625) q[3];
sx q[3];
rz(3.1282792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3964316) q[2];
sx q[2];
rz(-2.4691984) q[2];
sx q[2];
rz(-2.3154837) q[2];
rz(-2.7769026) q[3];
sx q[3];
rz(-1.6148022) q[3];
sx q[3];
rz(-0.29432347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18225886) q[0];
sx q[0];
rz(-1.8987645) q[0];
sx q[0];
rz(-0.44573319) q[0];
rz(0.590473) q[1];
sx q[1];
rz(-1.0346552) q[1];
sx q[1];
rz(0.42764923) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2954461) q[0];
sx q[0];
rz(-1.6177735) q[0];
sx q[0];
rz(-2.0278005) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1600464) q[2];
sx q[2];
rz(-0.84651154) q[2];
sx q[2];
rz(1.2638448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97274894) q[1];
sx q[1];
rz(-0.66509897) q[1];
sx q[1];
rz(1.14023) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1366051) q[3];
sx q[3];
rz(-0.79646275) q[3];
sx q[3];
rz(2.4465268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6815971) q[2];
sx q[2];
rz(-0.57969105) q[2];
sx q[2];
rz(-1.0350636) q[2];
rz(-1.5486108) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(2.2523994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11397938) q[0];
sx q[0];
rz(-3.1359105) q[0];
sx q[0];
rz(0.11153829) q[0];
rz(1.5882209) q[1];
sx q[1];
rz(-0.40451834) q[1];
sx q[1];
rz(-2.9389971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3521393) q[0];
sx q[0];
rz(-1.7048262) q[0];
sx q[0];
rz(-1.171597) q[0];
rz(0.28891383) q[2];
sx q[2];
rz(-1.8560858) q[2];
sx q[2];
rz(-0.71277518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21661585) q[1];
sx q[1];
rz(-1.0362715) q[1];
sx q[1];
rz(2.288544) q[1];
rz(-pi) q[2];
rz(-0.59059322) q[3];
sx q[3];
rz(-1.9265129) q[3];
sx q[3];
rz(1.4416665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69278875) q[2];
sx q[2];
rz(-1.0705798) q[2];
sx q[2];
rz(-2.1602574) q[2];
rz(0.28686178) q[3];
sx q[3];
rz(-0.57932866) q[3];
sx q[3];
rz(0.19744344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9051179) q[0];
sx q[0];
rz(-0.050967) q[0];
sx q[0];
rz(-2.4531051) q[0];
rz(-0.14432898) q[1];
sx q[1];
rz(-1.9537484) q[1];
sx q[1];
rz(0.7483288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2331325) q[0];
sx q[0];
rz(-1.5262506) q[0];
sx q[0];
rz(2.4219803) q[0];
rz(-1.2690721) q[2];
sx q[2];
rz(-1.8278484) q[2];
sx q[2];
rz(1.014705) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0341558) q[1];
sx q[1];
rz(-1.7567051) q[1];
sx q[1];
rz(0.30349813) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2714628) q[3];
sx q[3];
rz(-1.8091473) q[3];
sx q[3];
rz(0.42422653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.051108483) q[2];
sx q[2];
rz(-1.8277617) q[2];
sx q[2];
rz(-1.2468106) q[2];
rz(3.0221853) q[3];
sx q[3];
rz(-1.2179008) q[3];
sx q[3];
rz(0.30521211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7684105) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(2.5266732) q[0];
rz(2.3070295) q[1];
sx q[1];
rz(-1.7465697) q[1];
sx q[1];
rz(-0.83647299) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30395384) q[0];
sx q[0];
rz(-1.5364354) q[0];
sx q[0];
rz(-0.0092577309) q[0];
rz(-pi) q[1];
rz(-0.70186285) q[2];
sx q[2];
rz(-0.80282513) q[2];
sx q[2];
rz(2.7868164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3611722) q[1];
sx q[1];
rz(-0.57406146) q[1];
sx q[1];
rz(-1.3657544) q[1];
rz(-0.74252601) q[3];
sx q[3];
rz(-2.1603185) q[3];
sx q[3];
rz(1.3878915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60928434) q[2];
sx q[2];
rz(-0.60695761) q[2];
sx q[2];
rz(-0.92030805) q[2];
rz(0.80095428) q[3];
sx q[3];
rz(-2.6194173) q[3];
sx q[3];
rz(2.8029158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9698708) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(2.875476) q[0];
rz(-1.4316106) q[1];
sx q[1];
rz(-1.1191198) q[1];
sx q[1];
rz(-2.5439579) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1640329) q[0];
sx q[0];
rz(-0.64326875) q[0];
sx q[0];
rz(-0.32001647) q[0];
rz(-2.9605537) q[2];
sx q[2];
rz(-1.5597876) q[2];
sx q[2];
rz(2.371108) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2955851) q[1];
sx q[1];
rz(-2.4425828) q[1];
sx q[1];
rz(3.0547804) q[1];
rz(-pi) q[2];
rz(-1.131874) q[3];
sx q[3];
rz(-2.7898129) q[3];
sx q[3];
rz(1.1949415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0133682) q[2];
sx q[2];
rz(-0.2672264) q[2];
sx q[2];
rz(0.5989778) q[2];
rz(0.61601764) q[3];
sx q[3];
rz(-2.7003761) q[3];
sx q[3];
rz(2.0986572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041158572) q[0];
sx q[0];
rz(-0.91223311) q[0];
sx q[0];
rz(0.31901932) q[0];
rz(3.0306385) q[1];
sx q[1];
rz(-1.3919421) q[1];
sx q[1];
rz(-0.14161938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6017766) q[0];
sx q[0];
rz(-2.4519073) q[0];
sx q[0];
rz(-2.9897571) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0194112) q[2];
sx q[2];
rz(-1.5670318) q[2];
sx q[2];
rz(-2.5906445) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95699043) q[1];
sx q[1];
rz(-2.5612368) q[1];
sx q[1];
rz(0.34837153) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0504248) q[3];
sx q[3];
rz(-1.4401575) q[3];
sx q[3];
rz(-1.4121233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.55613279) q[2];
sx q[2];
rz(-1.0254878) q[2];
sx q[2];
rz(0.11036135) q[2];
rz(0.50802556) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(1.9917816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13462774) q[0];
sx q[0];
rz(-0.65161324) q[0];
sx q[0];
rz(-1.2192669) q[0];
rz(0.02267516) q[1];
sx q[1];
rz(-2.2494648) q[1];
sx q[1];
rz(-2.5882744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1071725) q[0];
sx q[0];
rz(-1.4822472) q[0];
sx q[0];
rz(0.91219418) q[0];
rz(1.0108286) q[2];
sx q[2];
rz(-0.65245318) q[2];
sx q[2];
rz(0.7459695) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0670315) q[1];
sx q[1];
rz(-2.0843049) q[1];
sx q[1];
rz(1.1719955) q[1];
x q[2];
rz(-0.24711547) q[3];
sx q[3];
rz(-1.2809506) q[3];
sx q[3];
rz(2.4245976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40552178) q[2];
sx q[2];
rz(-2.5260479) q[2];
sx q[2];
rz(-1.7919398) q[2];
rz(-3.0810629) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(2.5133666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8796006) q[0];
sx q[0];
rz(-0.99402004) q[0];
sx q[0];
rz(-0.0066198786) q[0];
rz(0.63893843) q[1];
sx q[1];
rz(-0.21955755) q[1];
sx q[1];
rz(-0.42983291) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13794261) q[0];
sx q[0];
rz(-1.2009504) q[0];
sx q[0];
rz(1.1843286) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0581994) q[2];
sx q[2];
rz(-0.48658961) q[2];
sx q[2];
rz(1.270382) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.51825) q[1];
sx q[1];
rz(-2.3584705) q[1];
sx q[1];
rz(3.1108969) q[1];
x q[2];
rz(-1.5019526) q[3];
sx q[3];
rz(-0.9864583) q[3];
sx q[3];
rz(-2.7513954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2868353) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(-3.1268934) q[2];
rz(1.1356575) q[3];
sx q[3];
rz(-2.0292) q[3];
sx q[3];
rz(1.8358102) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056674615) q[0];
sx q[0];
rz(-0.25923964) q[0];
sx q[0];
rz(1.2925451) q[0];
rz(-2.9855285) q[1];
sx q[1];
rz(-0.81764692) q[1];
sx q[1];
rz(-0.83052105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6097858) q[0];
sx q[0];
rz(-1.4713108) q[0];
sx q[0];
rz(-2.2983944) q[0];
rz(1.8761329) q[2];
sx q[2];
rz(-0.80080253) q[2];
sx q[2];
rz(0.82894553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76354181) q[1];
sx q[1];
rz(-1.4908067) q[1];
sx q[1];
rz(1.3536118) q[1];
rz(-pi) q[2];
rz(0.36673185) q[3];
sx q[3];
rz(-1.7730797) q[3];
sx q[3];
rz(2.3679581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0326651) q[2];
sx q[2];
rz(-2.445161) q[2];
sx q[2];
rz(-0.44309524) q[2];
rz(2.9923934) q[3];
sx q[3];
rz(-0.34330338) q[3];
sx q[3];
rz(-2.8086015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14015848) q[0];
sx q[0];
rz(-2.5509111) q[0];
sx q[0];
rz(2.1668707) q[0];
rz(-0.99826605) q[1];
sx q[1];
rz(-1.2405735) q[1];
sx q[1];
rz(2.1539198) q[1];
rz(-0.34299359) q[2];
sx q[2];
rz(-1.7609114) q[2];
sx q[2];
rz(-1.5821725) q[2];
rz(-2.0988437) q[3];
sx q[3];
rz(-2.4059912) q[3];
sx q[3];
rz(2.2148377) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
