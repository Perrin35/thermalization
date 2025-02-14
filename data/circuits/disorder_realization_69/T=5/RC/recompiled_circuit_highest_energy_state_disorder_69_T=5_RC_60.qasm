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
rz(0.34613553) q[0];
sx q[0];
rz(-1.5032285) q[0];
sx q[0];
rz(0.7315973) q[0];
rz(0.44959623) q[1];
sx q[1];
rz(-2.9540588) q[1];
sx q[1];
rz(-2.7076758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.514115) q[0];
sx q[0];
rz(-1.5356931) q[0];
sx q[0];
rz(2.3112464) q[0];
rz(-0.5753168) q[2];
sx q[2];
rz(-2.0849209) q[2];
sx q[2];
rz(-2.7935087) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0632532) q[1];
sx q[1];
rz(-1.9615478) q[1];
sx q[1];
rz(2.2558252) q[1];
x q[2];
rz(1.7619261) q[3];
sx q[3];
rz(-1.6014631) q[3];
sx q[3];
rz(0.14740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.71243858) q[2];
sx q[2];
rz(-1.466208) q[2];
sx q[2];
rz(-1.0038556) q[2];
rz(-0.53576523) q[3];
sx q[3];
rz(-0.86447132) q[3];
sx q[3];
rz(3.1404449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4681604) q[0];
sx q[0];
rz(-2.4591481) q[0];
sx q[0];
rz(-2.4144507) q[0];
rz(-0.06761059) q[1];
sx q[1];
rz(-1.7865684) q[1];
sx q[1];
rz(-2.0179857) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7794117) q[0];
sx q[0];
rz(-1.0955278) q[0];
sx q[0];
rz(-1.8888461) q[0];
x q[1];
rz(-3.0762663) q[2];
sx q[2];
rz(-1.4484754) q[2];
sx q[2];
rz(0.1647915) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.935168) q[1];
sx q[1];
rz(-1.4713396) q[1];
sx q[1];
rz(1.2897435) q[1];
rz(-pi) q[2];
rz(2.0455849) q[3];
sx q[3];
rz(-1.1752306) q[3];
sx q[3];
rz(-2.6717348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2850538) q[2];
sx q[2];
rz(-1.5679789) q[2];
sx q[2];
rz(0.48173586) q[2];
rz(0.5528062) q[3];
sx q[3];
rz(-1.0779287) q[3];
sx q[3];
rz(3.0520181) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32473096) q[0];
sx q[0];
rz(-2.6664) q[0];
sx q[0];
rz(0.09356308) q[0];
rz(-1.7612673) q[1];
sx q[1];
rz(-0.9535791) q[1];
sx q[1];
rz(2.5043452) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4555511) q[0];
sx q[0];
rz(-1.0555869) q[0];
sx q[0];
rz(-2.0705219) q[0];
rz(-pi) q[1];
rz(2.8222476) q[2];
sx q[2];
rz(-2.508906) q[2];
sx q[2];
rz(1.8377395) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4658424) q[1];
sx q[1];
rz(-1.9681962) q[1];
sx q[1];
rz(3.1134362) q[1];
rz(-pi) q[2];
rz(2.8508857) q[3];
sx q[3];
rz(-1.6120496) q[3];
sx q[3];
rz(0.34063971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9351585) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(0.68378249) q[2];
rz(0.23009662) q[3];
sx q[3];
rz(-1.4068539) q[3];
sx q[3];
rz(-0.42207119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(3.0919331) q[0];
sx q[0];
rz(-2.6347418) q[0];
sx q[0];
rz(1.7012713) q[0];
rz(0.47538844) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(-2.5604274) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6027045) q[0];
sx q[0];
rz(-1.4853444) q[0];
sx q[0];
rz(0.84792015) q[0];
rz(2.0330795) q[2];
sx q[2];
rz(-1.7812742) q[2];
sx q[2];
rz(2.0049948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95422259) q[1];
sx q[1];
rz(-1.1130325) q[1];
sx q[1];
rz(-1.2332031) q[1];
rz(-pi) q[2];
rz(2.851296) q[3];
sx q[3];
rz(-2.6321342) q[3];
sx q[3];
rz(0.51028937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0541957) q[2];
sx q[2];
rz(-0.18614686) q[2];
sx q[2];
rz(2.6206214) q[2];
rz(1.4969131) q[3];
sx q[3];
rz(-1.4119586) q[3];
sx q[3];
rz(-1.514667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5533376) q[0];
sx q[0];
rz(-2.8659358) q[0];
sx q[0];
rz(1.1302554) q[0];
rz(-2.0626119) q[1];
sx q[1];
rz(-1.9422453) q[1];
sx q[1];
rz(-2.030453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1905907) q[0];
sx q[0];
rz(-1.4050806) q[0];
sx q[0];
rz(-1.1072876) q[0];
rz(-pi) q[1];
rz(2.6118641) q[2];
sx q[2];
rz(-1.6090983) q[2];
sx q[2];
rz(1.8977579) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20566733) q[1];
sx q[1];
rz(-1.7280518) q[1];
sx q[1];
rz(-2.2155511) q[1];
x q[2];
rz(2.2419315) q[3];
sx q[3];
rz(-2.3658345) q[3];
sx q[3];
rz(0.44246261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7378716) q[2];
sx q[2];
rz(-0.51509866) q[2];
sx q[2];
rz(2.1007288) q[2];
rz(-0.8199842) q[3];
sx q[3];
rz(-1.9085725) q[3];
sx q[3];
rz(-0.6238873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89989221) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(2.4851121) q[0];
rz(1.3735324) q[1];
sx q[1];
rz(-2.3479925) q[1];
sx q[1];
rz(2.0097282) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8478125) q[0];
sx q[0];
rz(-0.72350693) q[0];
sx q[0];
rz(2.840994) q[0];
x q[1];
rz(2.4651338) q[2];
sx q[2];
rz(-2.0148811) q[2];
sx q[2];
rz(1.9537587) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0814077) q[1];
sx q[1];
rz(-2.2814732) q[1];
sx q[1];
rz(-3.1050472) q[1];
rz(-0.15040654) q[3];
sx q[3];
rz(-1.5751221) q[3];
sx q[3];
rz(-1.9459917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.446283) q[2];
sx q[2];
rz(-2.535847) q[2];
sx q[2];
rz(2.821935) q[2];
rz(-2.9122635) q[3];
sx q[3];
rz(-2.1181483) q[3];
sx q[3];
rz(-2.2314609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0730154) q[0];
sx q[0];
rz(-1.9972766) q[0];
sx q[0];
rz(-1.9101494) q[0];
rz(0.07864174) q[1];
sx q[1];
rz(-1.6784724) q[1];
sx q[1];
rz(-2.9873649) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.70792) q[0];
sx q[0];
rz(-1.3929954) q[0];
sx q[0];
rz(2.7975797) q[0];
rz(-2.8381589) q[2];
sx q[2];
rz(-1.9643485) q[2];
sx q[2];
rz(-2.9224666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5104766) q[1];
sx q[1];
rz(-2.0492024) q[1];
sx q[1];
rz(-1.9892742) q[1];
rz(-pi) q[2];
rz(1.7835791) q[3];
sx q[3];
rz(-1.9733505) q[3];
sx q[3];
rz(-2.8496131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8057033) q[2];
sx q[2];
rz(-1.615639) q[2];
sx q[2];
rz(-1.6501144) q[2];
rz(-1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(1.5707877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0989646) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(1.6866823) q[0];
rz(-0.97575724) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(-1.8029433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2678626) q[0];
sx q[0];
rz(-0.81276816) q[0];
sx q[0];
rz(1.0229179) q[0];
rz(-1.0089802) q[2];
sx q[2];
rz(-2.564401) q[2];
sx q[2];
rz(1.4711589) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6852831) q[1];
sx q[1];
rz(-1.3585618) q[1];
sx q[1];
rz(0.93385277) q[1];
rz(-0.46316163) q[3];
sx q[3];
rz(-2.4235631) q[3];
sx q[3];
rz(-2.5754186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6300388) q[2];
sx q[2];
rz(-2.126882) q[2];
sx q[2];
rz(-2.3021728) q[2];
rz(1.9073585) q[3];
sx q[3];
rz(-1.0890361) q[3];
sx q[3];
rz(-0.031410005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79427528) q[0];
sx q[0];
rz(-1.7592156) q[0];
sx q[0];
rz(2.7224702) q[0];
rz(-1.6436815) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(2.7083414) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99892761) q[0];
sx q[0];
rz(-1.6892528) q[0];
sx q[0];
rz(2.0303557) q[0];
rz(2.4845759) q[2];
sx q[2];
rz(-1.581356) q[2];
sx q[2];
rz(-2.8274909) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8318787) q[1];
sx q[1];
rz(-1.0597178) q[1];
sx q[1];
rz(2.8161777) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3550536) q[3];
sx q[3];
rz(-2.1387055) q[3];
sx q[3];
rz(-1.9892094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7822632) q[2];
sx q[2];
rz(-1.9908345) q[2];
sx q[2];
rz(2.9456054) q[2];
rz(2.0187142) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(-1.6897197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48542431) q[0];
sx q[0];
rz(-2.4805785) q[0];
sx q[0];
rz(2.0458903) q[0];
rz(-2.6307259) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(2.9313472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4097262) q[0];
sx q[0];
rz(-1.3748226) q[0];
sx q[0];
rz(-0.98742296) q[0];
x q[1];
rz(0.11504193) q[2];
sx q[2];
rz(-2.2224073) q[2];
sx q[2];
rz(0.63331214) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0970891) q[1];
sx q[1];
rz(-2.0538446) q[1];
sx q[1];
rz(-2.0562857) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5417762) q[3];
sx q[3];
rz(-1.2251523) q[3];
sx q[3];
rz(1.7315854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.11253396) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(1.0851592) q[2];
rz(2.8339913) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(-2.4881081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6899684) q[0];
sx q[0];
rz(-1.987048) q[0];
sx q[0];
rz(-1.2183627) q[0];
rz(2.4901509) q[1];
sx q[1];
rz(-0.94964288) q[1];
sx q[1];
rz(-2.3655187) q[1];
rz(0.015301558) q[2];
sx q[2];
rz(-2.2616784) q[2];
sx q[2];
rz(-0.26605284) q[2];
rz(-0.93919803) q[3];
sx q[3];
rz(-1.152335) q[3];
sx q[3];
rz(-1.4622968) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
