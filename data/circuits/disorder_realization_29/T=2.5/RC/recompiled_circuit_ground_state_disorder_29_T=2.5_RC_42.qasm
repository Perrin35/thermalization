OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.27958265) q[0];
sx q[0];
rz(-2.6065338) q[0];
sx q[0];
rz(-0.23393272) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(4.5555727) q[1];
sx q[1];
rz(11.04784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5748225) q[0];
sx q[0];
rz(-1.260551) q[0];
sx q[0];
rz(2.6022018) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1359623) q[2];
sx q[2];
rz(-1.7644492) q[2];
sx q[2];
rz(2.6676712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9649992) q[1];
sx q[1];
rz(-2.1273344) q[1];
sx q[1];
rz(2.2189369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9085479) q[3];
sx q[3];
rz(-1.9534142) q[3];
sx q[3];
rz(-2.1688714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5412377) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(-0.8055299) q[2];
rz(2.7496998) q[3];
sx q[3];
rz(-1.7854179) q[3];
sx q[3];
rz(-0.75683561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517752) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(0.43847325) q[0];
rz(-1.8665727) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(3.0335887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1718801) q[0];
sx q[0];
rz(-1.6816655) q[0];
sx q[0];
rz(-3.1165439) q[0];
rz(0.82307016) q[2];
sx q[2];
rz(-1.0727173) q[2];
sx q[2];
rz(-0.66253875) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0569563) q[1];
sx q[1];
rz(-2.8743187) q[1];
sx q[1];
rz(-1.8932503) q[1];
x q[2];
rz(-1.2233954) q[3];
sx q[3];
rz(-0.51651556) q[3];
sx q[3];
rz(1.4054738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7763623) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(-0.11621172) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(-2.3381086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47495833) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(-0.83589244) q[0];
rz(1.5180786) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(-1.2380884) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1956714) q[0];
sx q[0];
rz(-2.4950881) q[0];
sx q[0];
rz(-0.92574985) q[0];
rz(-pi) q[1];
rz(-2.798271) q[2];
sx q[2];
rz(-1.803071) q[2];
sx q[2];
rz(-2.6417123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9596531) q[1];
sx q[1];
rz(-1.5132071) q[1];
sx q[1];
rz(-1.4532386) q[1];
x q[2];
rz(-1.2830623) q[3];
sx q[3];
rz(-1.2997174) q[3];
sx q[3];
rz(-2.6379741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8572924) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(2.4227552) q[2];
rz(2.5607064) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(0.7287997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.523664) q[0];
sx q[0];
rz(-1.5939465) q[0];
sx q[0];
rz(2.0148328) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.6512066) q[1];
sx q[1];
rz(-1.0430956) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27985036) q[0];
sx q[0];
rz(-1.6011097) q[0];
sx q[0];
rz(1.3002731) q[0];
rz(-pi) q[1];
rz(2.1073694) q[2];
sx q[2];
rz(-1.3541647) q[2];
sx q[2];
rz(-2.72675) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.14341893) q[1];
sx q[1];
rz(-0.34356782) q[1];
sx q[1];
rz(1.0982537) q[1];
rz(-0.11123379) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(-1.8452132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2780693) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(-2.1045904) q[2];
rz(-0.95748025) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(-0.83468848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1400414) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(0.087604372) q[0];
rz(2.7032848) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(-1.3137438) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4632516) q[0];
sx q[0];
rz(-1.3366564) q[0];
sx q[0];
rz(-1.0749504) q[0];
rz(-pi) q[1];
rz(-1.472166) q[2];
sx q[2];
rz(-1.5336138) q[2];
sx q[2];
rz(-1.94869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6586761) q[1];
sx q[1];
rz(-0.19397846) q[1];
sx q[1];
rz(2.0493755) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99788061) q[3];
sx q[3];
rz(-0.30227236) q[3];
sx q[3];
rz(0.59461601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4801415) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(-2.5782149) q[2];
rz(-1.9208469) q[3];
sx q[3];
rz(-1.3910339) q[3];
sx q[3];
rz(-2.5105072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0193943) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(-0.011938183) q[0];
rz(2.7519233) q[1];
sx q[1];
rz(-2.4395112) q[1];
sx q[1];
rz(0.24519244) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58701506) q[0];
sx q[0];
rz(-2.3759807) q[0];
sx q[0];
rz(-0.9258201) q[0];
x q[1];
rz(-1.7952314) q[2];
sx q[2];
rz(-1.8230453) q[2];
sx q[2];
rz(-2.4750965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4575544) q[1];
sx q[1];
rz(-1.1361309) q[1];
sx q[1];
rz(-0.60575374) q[1];
rz(2.7390476) q[3];
sx q[3];
rz(-1.3005101) q[3];
sx q[3];
rz(-2.9077934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7616854) q[2];
sx q[2];
rz(-1.8319538) q[2];
sx q[2];
rz(-1.3690108) q[2];
rz(0.53064972) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(2.0556889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1928007) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(-2.8444994) q[0];
rz(-1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(-0.43103257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1566136) q[0];
sx q[0];
rz(-0.3604069) q[0];
sx q[0];
rz(-2.5094338) q[0];
rz(0.22054976) q[2];
sx q[2];
rz(-2.9052581) q[2];
sx q[2];
rz(2.5823809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.7073133) q[1];
sx q[1];
rz(-0.61798862) q[1];
sx q[1];
rz(0.98843482) q[1];
x q[2];
rz(0.2934612) q[3];
sx q[3];
rz(-1.0612504) q[3];
sx q[3];
rz(-2.5603385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83583528) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(0.26710278) q[2];
rz(0.032020656) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(-0.10572461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778462) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(-2.5015976) q[0];
rz(0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(-1.4124426) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0873358) q[0];
sx q[0];
rz(-2.1001513) q[0];
sx q[0];
rz(2.0482778) q[0];
x q[1];
rz(-1.6663297) q[2];
sx q[2];
rz(-0.49411202) q[2];
sx q[2];
rz(-1.4612559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0113653) q[1];
sx q[1];
rz(-1.2308569) q[1];
sx q[1];
rz(-1.4814427) q[1];
rz(-pi) q[2];
rz(-1.1453923) q[3];
sx q[3];
rz(-1.2266876) q[3];
sx q[3];
rz(-2.1895125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5414446) q[2];
sx q[2];
rz(-1.7057799) q[2];
sx q[2];
rz(-2.7195462) q[2];
rz(-0.47354928) q[3];
sx q[3];
rz(-0.62255064) q[3];
sx q[3];
rz(2.9464909) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710829) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(1.7899845) q[0];
rz(2.8260258) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(-1.9884761) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89698541) q[0];
sx q[0];
rz(-3.0685406) q[0];
sx q[0];
rz(2.3417419) q[0];
x q[1];
rz(0.53535992) q[2];
sx q[2];
rz(-1.6125814) q[2];
sx q[2];
rz(-0.9524065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7866384) q[1];
sx q[1];
rz(-1.4435205) q[1];
sx q[1];
rz(-0.58357012) q[1];
rz(-pi) q[2];
rz(2.7615776) q[3];
sx q[3];
rz(-1.9560062) q[3];
sx q[3];
rz(0.64982254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6844668) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(0.37377629) q[2];
rz(-1.9155546) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(1.3396243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(-1.8208338) q[0];
rz(-2.2044115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(-0.46930596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3166817) q[0];
sx q[0];
rz(-1.4535722) q[0];
sx q[0];
rz(0.076923142) q[0];
x q[1];
rz(-1.2861757) q[2];
sx q[2];
rz(-2.5409343) q[2];
sx q[2];
rz(2.9702368) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.093897029) q[1];
sx q[1];
rz(-2.9710431) q[1];
sx q[1];
rz(-0.49255126) q[1];
x q[2];
rz(-2.8402249) q[3];
sx q[3];
rz(-1.3346938) q[3];
sx q[3];
rz(2.2599041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70539537) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(-1.9311284) q[2];
rz(2.5734899) q[3];
sx q[3];
rz(-1.7567239) q[3];
sx q[3];
rz(2.1657522) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217011) q[0];
sx q[0];
rz(-2.2911063) q[0];
sx q[0];
rz(-1.6794857) q[0];
rz(-0.89314356) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(0.90441119) q[2];
sx q[2];
rz(-2.7795962) q[2];
sx q[2];
rz(-2.4301651) q[2];
rz(-2.1956069) q[3];
sx q[3];
rz(-2.47249) q[3];
sx q[3];
rz(-2.5853018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
