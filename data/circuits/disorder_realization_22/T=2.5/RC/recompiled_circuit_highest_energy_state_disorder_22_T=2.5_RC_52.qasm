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
rz(0.15777388) q[0];
sx q[0];
rz(-1.1717492) q[0];
sx q[0];
rz(-1.8921312) q[0];
rz(3.0013822) q[1];
sx q[1];
rz(-1.4445211) q[1];
sx q[1];
rz(-0.061847774) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88647283) q[0];
sx q[0];
rz(-0.97772163) q[0];
sx q[0];
rz(0.61086872) q[0];
rz(-pi) q[1];
rz(3.1377162) q[2];
sx q[2];
rz(-0.090423294) q[2];
sx q[2];
rz(-3.0341482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7951491) q[1];
sx q[1];
rz(-0.7278053) q[1];
sx q[1];
rz(-1.7176571) q[1];
x q[2];
rz(-2.0560045) q[3];
sx q[3];
rz(-0.74923979) q[3];
sx q[3];
rz(0.10252122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7141815) q[2];
sx q[2];
rz(-0.84250557) q[2];
sx q[2];
rz(-0.25644914) q[2];
rz(0.17476684) q[3];
sx q[3];
rz(-1.1656961) q[3];
sx q[3];
rz(2.3037361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3278219) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-0.12369618) q[0];
rz(0.23873121) q[1];
sx q[1];
rz(-0.78014603) q[1];
sx q[1];
rz(-3.0366268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5803842) q[0];
sx q[0];
rz(-0.7155025) q[0];
sx q[0];
rz(1.3405878) q[0];
rz(-pi) q[1];
rz(1.8027854) q[2];
sx q[2];
rz(-1.9769443) q[2];
sx q[2];
rz(-2.7628203) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8387977) q[1];
sx q[1];
rz(-1.0775591) q[1];
sx q[1];
rz(2.8431945) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4363748) q[3];
sx q[3];
rz(-1.2949484) q[3];
sx q[3];
rz(2.5666756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3671942) q[2];
sx q[2];
rz(-1.4906733) q[2];
sx q[2];
rz(-1.4614089) q[2];
rz(-1.8558308) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(-2.8694966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1478145) q[0];
sx q[0];
rz(-0.21037924) q[0];
sx q[0];
rz(-1.0149957) q[0];
rz(2.2442832) q[1];
sx q[1];
rz(-1.4941447) q[1];
sx q[1];
rz(-2.4651333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2059442) q[0];
sx q[0];
rz(-0.95399154) q[0];
sx q[0];
rz(-0.031868462) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6482267) q[2];
sx q[2];
rz(-1.0037494) q[2];
sx q[2];
rz(2.9285448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3520917) q[1];
sx q[1];
rz(-1.3695696) q[1];
sx q[1];
rz(-1.9580864) q[1];
rz(1.5900988) q[3];
sx q[3];
rz(-2.3870644) q[3];
sx q[3];
rz(1.7280662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93495381) q[2];
sx q[2];
rz(-0.97524869) q[2];
sx q[2];
rz(0.76033956) q[2];
rz(1.2065411) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(-0.84347239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43561414) q[0];
sx q[0];
rz(-2.5137081) q[0];
sx q[0];
rz(-1.5950369) q[0];
rz(0.71890038) q[1];
sx q[1];
rz(-1.8214106) q[1];
sx q[1];
rz(-2.9100606) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084615413) q[0];
sx q[0];
rz(-1.1418793) q[0];
sx q[0];
rz(-0.57022988) q[0];
rz(-pi) q[1];
rz(0.79265742) q[2];
sx q[2];
rz(-1.5680299) q[2];
sx q[2];
rz(-2.1994928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0147103) q[1];
sx q[1];
rz(-1.0877123) q[1];
sx q[1];
rz(-1.4654069) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9636964) q[3];
sx q[3];
rz(-2.1390954) q[3];
sx q[3];
rz(1.5400122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4270758) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(-2.5330949) q[2];
rz(0.19615873) q[3];
sx q[3];
rz(-1.4213296) q[3];
sx q[3];
rz(-0.99854809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2732368) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(-0.77907816) q[0];
rz(-2.211606) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(1.1735865) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54544696) q[0];
sx q[0];
rz(-1.9088424) q[0];
sx q[0];
rz(3.1374187) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0069238) q[2];
sx q[2];
rz(-2.1818518) q[2];
sx q[2];
rz(-0.5765644) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8820937) q[1];
sx q[1];
rz(-0.32392392) q[1];
sx q[1];
rz(-1.7739576) q[1];
rz(-pi) q[2];
rz(2.7125732) q[3];
sx q[3];
rz(-1.7439505) q[3];
sx q[3];
rz(2.8575983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.32213) q[2];
sx q[2];
rz(-0.16699114) q[2];
sx q[2];
rz(-0.67506153) q[2];
rz(-0.86554646) q[3];
sx q[3];
rz(-1.6277438) q[3];
sx q[3];
rz(1.2077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38782564) q[0];
sx q[0];
rz(-3.1238811) q[0];
sx q[0];
rz(-2.2702763) q[0];
rz(2.9246092) q[1];
sx q[1];
rz(-1.5996409) q[1];
sx q[1];
rz(-1.8035696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34505358) q[0];
sx q[0];
rz(-0.65700475) q[0];
sx q[0];
rz(-3.0057602) q[0];
rz(-pi) q[1];
rz(-1.7291405) q[2];
sx q[2];
rz(-0.73190231) q[2];
sx q[2];
rz(-0.5754234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2830401) q[1];
sx q[1];
rz(-1.6190908) q[1];
sx q[1];
rz(-1.4273391) q[1];
rz(0.48465259) q[3];
sx q[3];
rz(-1.402087) q[3];
sx q[3];
rz(1.3107306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.013082144) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(-2.3415372) q[2];
rz(2.1498146) q[3];
sx q[3];
rz(-0.9477152) q[3];
sx q[3];
rz(0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8868788) q[0];
sx q[0];
rz(-0.44097057) q[0];
sx q[0];
rz(0.15175858) q[0];
rz(-1.4166547) q[1];
sx q[1];
rz(-0.85985008) q[1];
sx q[1];
rz(-0.83622611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72513103) q[0];
sx q[0];
rz(-0.90354474) q[0];
sx q[0];
rz(-0.50531549) q[0];
rz(-pi) q[1];
rz(0.57668893) q[2];
sx q[2];
rz(-1.5637914) q[2];
sx q[2];
rz(-3.0617065) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9735896) q[1];
sx q[1];
rz(-2.8781946) q[1];
sx q[1];
rz(1.2348639) q[1];
rz(-pi) q[2];
rz(-2.8854675) q[3];
sx q[3];
rz(-0.47893347) q[3];
sx q[3];
rz(-2.0253945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4862711) q[2];
sx q[2];
rz(-3.0781015) q[2];
sx q[2];
rz(3.0625694) q[2];
rz(0.71845734) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(0.92459905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35543168) q[0];
sx q[0];
rz(-0.59816718) q[0];
sx q[0];
rz(0.86791903) q[0];
rz(2.4527841) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(-2.205663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0999231) q[0];
sx q[0];
rz(-1.7522305) q[0];
sx q[0];
rz(0.52930852) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2555588) q[2];
sx q[2];
rz(-0.5109957) q[2];
sx q[2];
rz(-3.0877132) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1327428) q[1];
sx q[1];
rz(-0.32948438) q[1];
sx q[1];
rz(1.7368172) q[1];
rz(3.1400263) q[3];
sx q[3];
rz(-0.74347444) q[3];
sx q[3];
rz(1.8494128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1615289) q[2];
sx q[2];
rz(-2.4825725) q[2];
sx q[2];
rz(0.25679437) q[2];
rz(-2.1234546) q[3];
sx q[3];
rz(-1.1192106) q[3];
sx q[3];
rz(0.97091466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37860206) q[0];
sx q[0];
rz(-0.55736962) q[0];
sx q[0];
rz(-2.2879404) q[0];
rz(1.8866106) q[1];
sx q[1];
rz(-1.1458784) q[1];
sx q[1];
rz(2.8010211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69215323) q[0];
sx q[0];
rz(-1.3763577) q[0];
sx q[0];
rz(-2.9783335) q[0];
x q[1];
rz(-2.4677327) q[2];
sx q[2];
rz(-1.7654716) q[2];
sx q[2];
rz(-1.0985634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53696991) q[1];
sx q[1];
rz(-0.83965644) q[1];
sx q[1];
rz(-0.62418681) q[1];
x q[2];
rz(-0.27517516) q[3];
sx q[3];
rz(-2.0180297) q[3];
sx q[3];
rz(2.3768611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42369947) q[2];
sx q[2];
rz(-1.3807978) q[2];
sx q[2];
rz(-2.763486) q[2];
rz(-0.2374436) q[3];
sx q[3];
rz(-1.0708555) q[3];
sx q[3];
rz(2.8606991) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045192748) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(-2.4943446) q[0];
rz(-1.1680394) q[1];
sx q[1];
rz(-2.2868575) q[1];
sx q[1];
rz(2.7580269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3824359) q[0];
sx q[0];
rz(-1.6014227) q[0];
sx q[0];
rz(1.5492065) q[0];
x q[1];
rz(-2.6099989) q[2];
sx q[2];
rz(-1.8583014) q[2];
sx q[2];
rz(-0.93898857) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4870461) q[1];
sx q[1];
rz(-0.33815835) q[1];
sx q[1];
rz(-0.66429331) q[1];
rz(3.0342919) q[3];
sx q[3];
rz(-2.0837874) q[3];
sx q[3];
rz(-1.9657621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26455227) q[2];
sx q[2];
rz(-2.2525747) q[2];
sx q[2];
rz(-1.5569347) q[2];
rz(1.6282188) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(-0.42678601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601111) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(-1.9602641) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(-1.0931821) q[2];
sx q[2];
rz(-1.7872767) q[2];
sx q[2];
rz(1.7959303) q[2];
rz(-0.48458002) q[3];
sx q[3];
rz(-1.598834) q[3];
sx q[3];
rz(-1.9091189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
