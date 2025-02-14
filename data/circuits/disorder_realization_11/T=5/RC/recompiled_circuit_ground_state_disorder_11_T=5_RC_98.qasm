OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.017555822) q[0];
sx q[0];
rz(3.2942009) q[0];
sx q[0];
rz(9.6080544) q[0];
rz(-2.2493989) q[1];
sx q[1];
rz(-2.0018938) q[1];
sx q[1];
rz(-2.5633864) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728672) q[0];
sx q[0];
rz(-1.5807165) q[0];
sx q[0];
rz(-3.0501151) q[0];
rz(-pi) q[1];
rz(-0.033138795) q[2];
sx q[2];
rz(-1.4862747) q[2];
sx q[2];
rz(-3.0042841) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6975054) q[1];
sx q[1];
rz(-2.0365391) q[1];
sx q[1];
rz(2.5986141) q[1];
rz(-pi) q[2];
rz(1.4981323) q[3];
sx q[3];
rz(-2.6748219) q[3];
sx q[3];
rz(2.6456986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1766498) q[2];
sx q[2];
rz(-0.86805934) q[2];
sx q[2];
rz(-3.0435666) q[2];
rz(2.3227504) q[3];
sx q[3];
rz(-0.10418532) q[3];
sx q[3];
rz(-2.5067743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1408511) q[0];
sx q[0];
rz(-2.8241557) q[0];
sx q[0];
rz(2.4733518) q[0];
rz(2.5375598) q[1];
sx q[1];
rz(-3.1112473) q[1];
sx q[1];
rz(-2.277453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1304431) q[0];
sx q[0];
rz(-0.98447039) q[0];
sx q[0];
rz(-1.4545813) q[0];
rz(-2.2496944) q[2];
sx q[2];
rz(-1.9111655) q[2];
sx q[2];
rz(0.75770411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27990231) q[1];
sx q[1];
rz(-0.30843014) q[1];
sx q[1];
rz(2.4753086) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92738028) q[3];
sx q[3];
rz(-2.815554) q[3];
sx q[3];
rz(0.044635208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16022564) q[2];
sx q[2];
rz(-1.9619433) q[2];
sx q[2];
rz(-0.78440624) q[2];
rz(-1.9085599) q[3];
sx q[3];
rz(-2.8630856) q[3];
sx q[3];
rz(1.9306785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8104372) q[0];
sx q[0];
rz(-1.5758608) q[0];
sx q[0];
rz(-1.021215) q[0];
rz(-1.637623) q[1];
sx q[1];
rz(-2.8228788) q[1];
sx q[1];
rz(2.5655897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2191781) q[0];
sx q[0];
rz(-1.4542581) q[0];
sx q[0];
rz(-0.4187421) q[0];
rz(-0.25635135) q[2];
sx q[2];
rz(-1.1927989) q[2];
sx q[2];
rz(-2.5963654) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6627548) q[1];
sx q[1];
rz(-0.75482268) q[1];
sx q[1];
rz(-2.635582) q[1];
rz(-pi) q[2];
rz(-0.66815362) q[3];
sx q[3];
rz(-0.34121554) q[3];
sx q[3];
rz(1.4618272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3699428) q[2];
sx q[2];
rz(-0.45054951) q[2];
sx q[2];
rz(0.096262781) q[2];
rz(0.43109584) q[3];
sx q[3];
rz(-0.96612203) q[3];
sx q[3];
rz(2.9935484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3189321) q[0];
sx q[0];
rz(-3.0156101) q[0];
sx q[0];
rz(-0.2952964) q[0];
rz(-2.2604306) q[1];
sx q[1];
rz(-0.22717871) q[1];
sx q[1];
rz(2.6491902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36400255) q[0];
sx q[0];
rz(-0.26960281) q[0];
sx q[0];
rz(2.2436938) q[0];
rz(-pi) q[1];
rz(3.0922331) q[2];
sx q[2];
rz(-1.4774895) q[2];
sx q[2];
rz(-2.9010584) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0379115) q[1];
sx q[1];
rz(-1.7010444) q[1];
sx q[1];
rz(1.9027846) q[1];
rz(-pi) q[2];
rz(1.9183228) q[3];
sx q[3];
rz(-1.9869958) q[3];
sx q[3];
rz(0.89194854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1262576) q[2];
sx q[2];
rz(-2.8110562) q[2];
sx q[2];
rz(-2.7988561) q[2];
rz(-0.98894173) q[3];
sx q[3];
rz(-1.9804695) q[3];
sx q[3];
rz(2.4435918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2547176) q[0];
sx q[0];
rz(-2.8509792) q[0];
sx q[0];
rz(-1.2683723) q[0];
rz(-0.36240029) q[1];
sx q[1];
rz(-0.86530322) q[1];
sx q[1];
rz(1.5579582) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9824955) q[0];
sx q[0];
rz(-1.0282059) q[0];
sx q[0];
rz(-2.1033903) q[0];
x q[1];
rz(2.1554016) q[2];
sx q[2];
rz(-0.3282983) q[2];
sx q[2];
rz(-2.3165348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6889188) q[1];
sx q[1];
rz(-2.5443379) q[1];
sx q[1];
rz(-0.026845499) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51934989) q[3];
sx q[3];
rz(-1.2548813) q[3];
sx q[3];
rz(-2.702932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6543988) q[2];
sx q[2];
rz(-2.0687658) q[2];
sx q[2];
rz(-2.2844592) q[2];
rz(-1.0355787) q[3];
sx q[3];
rz(-2.3643957) q[3];
sx q[3];
rz(-2.5100759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71984464) q[0];
sx q[0];
rz(-0.48374614) q[0];
sx q[0];
rz(2.8068722) q[0];
rz(2.1597629) q[1];
sx q[1];
rz(-1.5900759) q[1];
sx q[1];
rz(2.2614711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28583254) q[0];
sx q[0];
rz(-1.8632392) q[0];
sx q[0];
rz(0.16181914) q[0];
x q[1];
rz(0.99079972) q[2];
sx q[2];
rz(-0.31337619) q[2];
sx q[2];
rz(-0.51381451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.0047061027) q[1];
sx q[1];
rz(-2.0396173) q[1];
sx q[1];
rz(-0.15699082) q[1];
rz(-pi) q[2];
rz(0.91428925) q[3];
sx q[3];
rz(-0.56826545) q[3];
sx q[3];
rz(1.9075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52241391) q[2];
sx q[2];
rz(-2.652707) q[2];
sx q[2];
rz(1.4948814) q[2];
rz(1.304168) q[3];
sx q[3];
rz(-2.2924278) q[3];
sx q[3];
rz(-2.9554534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78231597) q[0];
sx q[0];
rz(-0.19141153) q[0];
sx q[0];
rz(-0.070847832) q[0];
rz(2.518636) q[1];
sx q[1];
rz(-2.7406335) q[1];
sx q[1];
rz(-0.43359044) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48468785) q[0];
sx q[0];
rz(-1.6015717) q[0];
sx q[0];
rz(-2.8367443) q[0];
rz(-1.5135399) q[2];
sx q[2];
rz(-1.861915) q[2];
sx q[2];
rz(1.2886485) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2411464) q[1];
sx q[1];
rz(-1.6668607) q[1];
sx q[1];
rz(0.32498951) q[1];
x q[2];
rz(2.904312) q[3];
sx q[3];
rz(-2.6868636) q[3];
sx q[3];
rz(0.90516289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5101461) q[2];
sx q[2];
rz(-2.6084709) q[2];
sx q[2];
rz(-0.637429) q[2];
rz(-2.0006477) q[3];
sx q[3];
rz(-0.44064042) q[3];
sx q[3];
rz(-2.2250037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7618074) q[0];
sx q[0];
rz(-0.26476911) q[0];
sx q[0];
rz(2.8763212) q[0];
rz(-0.028566407) q[1];
sx q[1];
rz(-0.18762372) q[1];
sx q[1];
rz(-2.2144337) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1409455) q[0];
sx q[0];
rz(-0.22613284) q[0];
sx q[0];
rz(-0.12238149) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4248542) q[2];
sx q[2];
rz(-1.2044014) q[2];
sx q[2];
rz(1.4595569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68009752) q[1];
sx q[1];
rz(-2.3465183) q[1];
sx q[1];
rz(2.9976588) q[1];
rz(-pi) q[2];
rz(0.42334564) q[3];
sx q[3];
rz(-1.6739427) q[3];
sx q[3];
rz(2.1684017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.022543) q[2];
sx q[2];
rz(-1.0485342) q[2];
sx q[2];
rz(2.0691464) q[2];
rz(0.33506814) q[3];
sx q[3];
rz(-3.0266893) q[3];
sx q[3];
rz(-0.2255628) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2304614) q[0];
sx q[0];
rz(-1.1536396) q[0];
sx q[0];
rz(-2.352584) q[0];
rz(2.543653) q[1];
sx q[1];
rz(-2.0409248) q[1];
sx q[1];
rz(2.423563) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5563468) q[0];
sx q[0];
rz(-1.5296773) q[0];
sx q[0];
rz(-1.7264191) q[0];
rz(-pi) q[1];
rz(-1.687019) q[2];
sx q[2];
rz(-1.5614484) q[2];
sx q[2];
rz(-2.7264894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2749697) q[1];
sx q[1];
rz(-2.1964745) q[1];
sx q[1];
rz(0.83557202) q[1];
x q[2];
rz(-0.17810693) q[3];
sx q[3];
rz(-1.6619887) q[3];
sx q[3];
rz(-0.9739463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52704063) q[2];
sx q[2];
rz(-0.25200945) q[2];
sx q[2];
rz(-1.7238114) q[2];
rz(-2.0333911) q[3];
sx q[3];
rz(-0.36644822) q[3];
sx q[3];
rz(0.38780701) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.336816) q[0];
sx q[0];
rz(-0.62938654) q[0];
sx q[0];
rz(0.25548536) q[0];
rz(-0.77087036) q[1];
sx q[1];
rz(-1.5522542) q[1];
sx q[1];
rz(1.5482056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436018) q[0];
sx q[0];
rz(-1.457347) q[0];
sx q[0];
rz(1.1313637) q[0];
rz(-pi) q[1];
rz(1.3012655) q[2];
sx q[2];
rz(-2.4204006) q[2];
sx q[2];
rz(2.1047956) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3774958) q[1];
sx q[1];
rz(-1.2974655) q[1];
sx q[1];
rz(0.39876858) q[1];
rz(-pi) q[2];
rz(1.19243) q[3];
sx q[3];
rz(-2.1244266) q[3];
sx q[3];
rz(2.7146482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2737246) q[2];
sx q[2];
rz(-0.74213433) q[2];
sx q[2];
rz(1.8871657) q[2];
rz(-0.54002386) q[3];
sx q[3];
rz(-0.050914474) q[3];
sx q[3];
rz(2.36256) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9191606) q[0];
sx q[0];
rz(-1.5615015) q[0];
sx q[0];
rz(-1.1175565) q[0];
rz(-0.22668214) q[1];
sx q[1];
rz(-0.10348987) q[1];
sx q[1];
rz(-1.6685974) q[1];
rz(-0.37353362) q[2];
sx q[2];
rz(-2.1185771) q[2];
sx q[2];
rz(2.6368903) q[2];
rz(2.7829783) q[3];
sx q[3];
rz(-2.3862541) q[3];
sx q[3];
rz(-2.42498) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
