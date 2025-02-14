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
rz(0.67233664) q[0];
sx q[0];
rz(5.0007102) q[0];
sx q[0];
rz(11.034427) q[0];
rz(1.2338282) q[1];
sx q[1];
rz(-1.9246074) q[1];
sx q[1];
rz(-1.4438862) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1506526) q[0];
sx q[0];
rz(-2.4661337) q[0];
sx q[0];
rz(-1.4399685) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9305761) q[2];
sx q[2];
rz(-1.8471247) q[2];
sx q[2];
rz(0.88283774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4593764) q[1];
sx q[1];
rz(-0.77118783) q[1];
sx q[1];
rz(2.5384063) q[1];
rz(3.1137556) q[3];
sx q[3];
rz(-1.2317766) q[3];
sx q[3];
rz(-2.8852685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3621138) q[2];
sx q[2];
rz(-1.9356091) q[2];
sx q[2];
rz(2.3625679) q[2];
rz(1.3075167) q[3];
sx q[3];
rz(-2.8179759) q[3];
sx q[3];
rz(1.3870846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55325145) q[0];
sx q[0];
rz(-1.7760176) q[0];
sx q[0];
rz(-2.8643082) q[0];
rz(-2.7514027) q[1];
sx q[1];
rz(-1.1017825) q[1];
sx q[1];
rz(-1.1276833) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4235315) q[0];
sx q[0];
rz(-3.09457) q[0];
sx q[0];
rz(-2.7014947) q[0];
rz(-pi) q[1];
rz(-1.4696944) q[2];
sx q[2];
rz(-1.8549457) q[2];
sx q[2];
rz(2.2997039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1192644) q[1];
sx q[1];
rz(-0.27909595) q[1];
sx q[1];
rz(-1.3420245) q[1];
rz(0.66172285) q[3];
sx q[3];
rz(-1.8764917) q[3];
sx q[3];
rz(2.2591285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4358501) q[2];
sx q[2];
rz(-0.59429344) q[2];
sx q[2];
rz(-1.1951813) q[2];
rz(-2.4863906) q[3];
sx q[3];
rz(-1.4846669) q[3];
sx q[3];
rz(0.5932194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6708267) q[0];
sx q[0];
rz(-2.4556181) q[0];
sx q[0];
rz(-0.92054787) q[0];
rz(1.8596733) q[1];
sx q[1];
rz(-2.0125407) q[1];
sx q[1];
rz(1.6887853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0272511) q[0];
sx q[0];
rz(-2.998861) q[0];
sx q[0];
rz(-2.4817075) q[0];
rz(-pi) q[1];
rz(-1.0017582) q[2];
sx q[2];
rz(-1.8382706) q[2];
sx q[2];
rz(-0.75050577) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3753759) q[1];
sx q[1];
rz(-1.2355243) q[1];
sx q[1];
rz(1.6250618) q[1];
x q[2];
rz(-1.6369197) q[3];
sx q[3];
rz(-1.6363314) q[3];
sx q[3];
rz(-1.6082482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.526223) q[2];
sx q[2];
rz(-1.0542032) q[2];
sx q[2];
rz(-2.7426381) q[2];
rz(0.51359549) q[3];
sx q[3];
rz(-0.31702888) q[3];
sx q[3];
rz(1.1494466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8154163) q[0];
sx q[0];
rz(-2.5626817) q[0];
sx q[0];
rz(1.2166566) q[0];
rz(-0.51219621) q[1];
sx q[1];
rz(-2.0307816) q[1];
sx q[1];
rz(-2.0373352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7889658) q[0];
sx q[0];
rz(-1.8471083) q[0];
sx q[0];
rz(-0.31310149) q[0];
rz(-pi) q[1];
rz(2.3338564) q[2];
sx q[2];
rz(-2.8931354) q[2];
sx q[2];
rz(-2.823644) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19266549) q[1];
sx q[1];
rz(-1.7671698) q[1];
sx q[1];
rz(3.0801303) q[1];
rz(-pi) q[2];
rz(-0.85569546) q[3];
sx q[3];
rz(-2.8036806) q[3];
sx q[3];
rz(-2.5355946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5451374) q[2];
sx q[2];
rz(-1.1918273) q[2];
sx q[2];
rz(-0.54984251) q[2];
rz(0.95692974) q[3];
sx q[3];
rz(-2.3423829) q[3];
sx q[3];
rz(-1.4889312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(1.4855708) q[0];
sx q[0];
rz(-0.85000426) q[0];
sx q[0];
rz(-2.3824084) q[0];
rz(-0.94866577) q[1];
sx q[1];
rz(-1.1439088) q[1];
sx q[1];
rz(2.8099828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0873957) q[0];
sx q[0];
rz(-0.83354649) q[0];
sx q[0];
rz(-1.1655318) q[0];
rz(-1.5344083) q[2];
sx q[2];
rz(-2.0452283) q[2];
sx q[2];
rz(-0.74288054) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79072551) q[1];
sx q[1];
rz(-1.7247883) q[1];
sx q[1];
rz(1.3920946) q[1];
rz(-pi) q[2];
rz(1.4344352) q[3];
sx q[3];
rz(-1.9788187) q[3];
sx q[3];
rz(1.3566164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18629508) q[2];
sx q[2];
rz(-2.1897327) q[2];
sx q[2];
rz(-2.1985506) q[2];
rz(-0.17213639) q[3];
sx q[3];
rz(-2.19682) q[3];
sx q[3];
rz(-0.22323639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73582369) q[0];
sx q[0];
rz(-2.3961234) q[0];
sx q[0];
rz(-1.5752342) q[0];
rz(-0.15815059) q[1];
sx q[1];
rz(-2.3731396) q[1];
sx q[1];
rz(2.2134773) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52735746) q[0];
sx q[0];
rz(-1.4422461) q[0];
sx q[0];
rz(0.30621333) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0177703) q[2];
sx q[2];
rz(-1.9352416) q[2];
sx q[2];
rz(-1.5096957) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0911407) q[1];
sx q[1];
rz(-0.71040043) q[1];
sx q[1];
rz(1.1732167) q[1];
rz(-0.32676231) q[3];
sx q[3];
rz(-2.770677) q[3];
sx q[3];
rz(-2.4075631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9485665) q[2];
sx q[2];
rz(-0.28634772) q[2];
sx q[2];
rz(-2.7890653) q[2];
rz(-2.7726717) q[3];
sx q[3];
rz(-2.5589294) q[3];
sx q[3];
rz(-2.2841456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.3567268) q[0];
sx q[0];
rz(-0.016473869) q[0];
sx q[0];
rz(-2.4705868) q[0];
rz(-3.0352133) q[1];
sx q[1];
rz(-2.2439067) q[1];
sx q[1];
rz(-2.5118714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9167921) q[0];
sx q[0];
rz(-0.19242254) q[0];
sx q[0];
rz(-3.0436361) q[0];
rz(-pi) q[1];
rz(0.81083576) q[2];
sx q[2];
rz(-0.5631643) q[2];
sx q[2];
rz(2.882618) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32618648) q[1];
sx q[1];
rz(-1.2478653) q[1];
sx q[1];
rz(1.6492483) q[1];
x q[2];
rz(1.2604146) q[3];
sx q[3];
rz(-2.2058626) q[3];
sx q[3];
rz(2.6657651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.070039198) q[2];
sx q[2];
rz(-0.89954251) q[2];
sx q[2];
rz(-2.1001935) q[2];
rz(-1.6085767) q[3];
sx q[3];
rz(-0.93419111) q[3];
sx q[3];
rz(1.9937203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94569412) q[0];
sx q[0];
rz(-2.4905289) q[0];
sx q[0];
rz(-0.4796637) q[0];
rz(3.116963) q[1];
sx q[1];
rz(-1.004091) q[1];
sx q[1];
rz(-3.0020795) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5845156) q[0];
sx q[0];
rz(-3.0026443) q[0];
sx q[0];
rz(0.90451972) q[0];
rz(-pi) q[1];
rz(-2.191338) q[2];
sx q[2];
rz(-0.59796158) q[2];
sx q[2];
rz(0.68956748) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3132726) q[1];
sx q[1];
rz(-1.7370151) q[1];
sx q[1];
rz(-0.66201026) q[1];
rz(-pi) q[2];
rz(1.516253) q[3];
sx q[3];
rz(-2.2404376) q[3];
sx q[3];
rz(3.0564133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8861683) q[2];
sx q[2];
rz(-1.5972127) q[2];
sx q[2];
rz(1.1523979) q[2];
rz(-1.7216916) q[3];
sx q[3];
rz(-2.4072188) q[3];
sx q[3];
rz(-1.3091807) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060147978) q[0];
sx q[0];
rz(-1.4652493) q[0];
sx q[0];
rz(-1.6434705) q[0];
rz(2.3807047) q[1];
sx q[1];
rz(-0.94327578) q[1];
sx q[1];
rz(1.6852185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13933548) q[0];
sx q[0];
rz(-0.59334785) q[0];
sx q[0];
rz(2.8171881) q[0];
rz(-pi) q[1];
rz(-1.782519) q[2];
sx q[2];
rz(-1.178732) q[2];
sx q[2];
rz(1.7535007) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.06157626) q[1];
sx q[1];
rz(-1.1760839) q[1];
sx q[1];
rz(1.030974) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1377148) q[3];
sx q[3];
rz(-1.9291546) q[3];
sx q[3];
rz(1.5556415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8196408) q[2];
sx q[2];
rz(-2.7072622) q[2];
sx q[2];
rz(-1.0158018) q[2];
rz(-2.2470233) q[3];
sx q[3];
rz(-1.8958478) q[3];
sx q[3];
rz(-2.3999124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1557409) q[0];
sx q[0];
rz(-2.4116801) q[0];
sx q[0];
rz(-0.99005449) q[0];
rz(0.92297018) q[1];
sx q[1];
rz(-1.762941) q[1];
sx q[1];
rz(-2.1754481) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21912757) q[0];
sx q[0];
rz(-0.12258633) q[0];
sx q[0];
rz(2.1845093) q[0];
x q[1];
rz(-0.43010148) q[2];
sx q[2];
rz(-1.9942339) q[2];
sx q[2];
rz(0.11030876) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6539972) q[1];
sx q[1];
rz(-1.8593809) q[1];
sx q[1];
rz(1.9777795) q[1];
x q[2];
rz(0.11364348) q[3];
sx q[3];
rz(-1.7688278) q[3];
sx q[3];
rz(-2.6012492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82603377) q[2];
sx q[2];
rz(-2.3580599) q[2];
sx q[2];
rz(-2.3890736) q[2];
rz(0.13713947) q[3];
sx q[3];
rz(-0.64118853) q[3];
sx q[3];
rz(-2.7916059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8990477) q[0];
sx q[0];
rz(-0.62472961) q[0];
sx q[0];
rz(-0.20150264) q[0];
rz(1.3626199) q[1];
sx q[1];
rz(-2.2687804) q[1];
sx q[1];
rz(-0.16558095) q[1];
rz(-1.0612442) q[2];
sx q[2];
rz(-0.10165086) q[2];
sx q[2];
rz(1.0063444) q[2];
rz(0.36214245) q[3];
sx q[3];
rz(-1.9531622) q[3];
sx q[3];
rz(1.2827557) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
