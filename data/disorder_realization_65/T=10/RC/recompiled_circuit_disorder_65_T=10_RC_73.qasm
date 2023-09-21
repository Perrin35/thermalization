OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(1.6488099) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(5.5881349) q[1];
sx q[1];
rz(10.499428) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2225392) q[0];
sx q[0];
rz(-0.018964501) q[0];
sx q[0];
rz(2.8595964) q[0];
rz(-pi) q[1];
rz(-1.4104112) q[2];
sx q[2];
rz(-2.3540902) q[2];
sx q[2];
rz(3.0410142) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5584471) q[1];
sx q[1];
rz(-1.6513283) q[1];
sx q[1];
rz(3.0196378) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69016506) q[3];
sx q[3];
rz(-0.65342045) q[3];
sx q[3];
rz(2.0131468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.82912123) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.6281698) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(-0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56869498) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(0.16076316) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(-1.6404023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06022913) q[0];
sx q[0];
rz(-1.4481359) q[0];
sx q[0];
rz(1.5457904) q[0];
rz(2.03962) q[2];
sx q[2];
rz(-1.8509682) q[2];
sx q[2];
rz(2.3902992) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0928287) q[1];
sx q[1];
rz(-1.7424912) q[1];
sx q[1];
rz(2.462639) q[1];
rz(1.5853911) q[3];
sx q[3];
rz(-2.7505034) q[3];
sx q[3];
rz(0.99816834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(1.1668011) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(-0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(1.5037781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30341879) q[0];
sx q[0];
rz(-1.5056599) q[0];
sx q[0];
rz(1.6584048) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90942803) q[2];
sx q[2];
rz(-0.69790188) q[2];
sx q[2];
rz(2.1991889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7392002) q[1];
sx q[1];
rz(-1.8715197) q[1];
sx q[1];
rz(1.6154358) q[1];
rz(2.2272439) q[3];
sx q[3];
rz(-1.4667061) q[3];
sx q[3];
rz(-0.57746938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(2.0898576) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.1490885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79608941) q[0];
sx q[0];
rz(-0.26608276) q[0];
sx q[0];
rz(-2.7059113) q[0];
rz(2.1637784) q[2];
sx q[2];
rz(-1.9359509) q[2];
sx q[2];
rz(2.1821373) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.541084) q[1];
sx q[1];
rz(-2.2693172) q[1];
sx q[1];
rz(-0.38703106) q[1];
rz(-pi) q[2];
rz(-1.9150919) q[3];
sx q[3];
rz(-0.85840423) q[3];
sx q[3];
rz(0.32574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2465308) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(2.4728298) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(1.0523798) q[0];
rz(1.2106238) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(2.2573684) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852172) q[0];
sx q[0];
rz(-1.3803122) q[0];
sx q[0];
rz(2.9604244) q[0];
rz(-2.7344633) q[2];
sx q[2];
rz(-1.1659157) q[2];
sx q[2];
rz(1.4332353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3430749) q[1];
sx q[1];
rz(-1.6371562) q[1];
sx q[1];
rz(2.9970478) q[1];
rz(-pi) q[2];
x q[2];
rz(2.103881) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(2.0562293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9050682) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(1.0926251) q[2];
rz(-0.24400273) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(2.3905335) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(2.172519) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5613865) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(-1.6798621) q[0];
rz(3.063383) q[2];
sx q[2];
rz(-0.71560301) q[2];
sx q[2];
rz(2.2323334) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8816996) q[1];
sx q[1];
rz(-2.0310146) q[1];
sx q[1];
rz(1.1020745) q[1];
rz(0.46521516) q[3];
sx q[3];
rz(-1.0930659) q[3];
sx q[3];
rz(2.5177237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(-0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(-1.4490022) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(-1.70599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3754417) q[0];
sx q[0];
rz(-1.0784082) q[0];
sx q[0];
rz(2.173461) q[0];
x q[1];
rz(0.018354015) q[2];
sx q[2];
rz(-1.5510501) q[2];
sx q[2];
rz(0.28723082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.095852921) q[1];
sx q[1];
rz(-0.47912712) q[1];
sx q[1];
rz(-2.2687001) q[1];
rz(0.88019754) q[3];
sx q[3];
rz(-0.62303632) q[3];
sx q[3];
rz(2.8043889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-3.0086349) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(-1.3652323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62194659) q[0];
sx q[0];
rz(-1.4953817) q[0];
sx q[0];
rz(-1.4343065) q[0];
x q[1];
rz(-1.6985967) q[2];
sx q[2];
rz(-1.5926077) q[2];
sx q[2];
rz(-2.4335361) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8147162) q[1];
sx q[1];
rz(-1.5183581) q[1];
sx q[1];
rz(0.53175064) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41994628) q[3];
sx q[3];
rz(-1.4010324) q[3];
sx q[3];
rz(-1.9672293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0662971) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-2.7189642) q[2];
rz(-2.1832809) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1388824) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(0.92393595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5101178) q[0];
sx q[0];
rz(-2.3259813) q[0];
sx q[0];
rz(-2.5328013) q[0];
x q[1];
rz(0.48859476) q[2];
sx q[2];
rz(-1.5813507) q[2];
sx q[2];
rz(0.34305629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82211557) q[1];
sx q[1];
rz(-1.9648106) q[1];
sx q[1];
rz(-2.134321) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1903673) q[3];
sx q[3];
rz(-2.1853672) q[3];
sx q[3];
rz(-2.8545477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.7473934) q[2];
rz(-1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(-2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37344638) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.753153) q[0];
rz(0.0069847981) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-3.1034234) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.105135) q[0];
sx q[0];
rz(-1.5441226) q[0];
sx q[0];
rz(-0.96203502) q[0];
rz(-pi) q[1];
rz(1.9534692) q[2];
sx q[2];
rz(-1.4854991) q[2];
sx q[2];
rz(1.485699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2529777) q[1];
sx q[1];
rz(-1.1974317) q[1];
sx q[1];
rz(0.044815973) q[1];
rz(0.84382236) q[3];
sx q[3];
rz(-2.3355964) q[3];
sx q[3];
rz(-2.3674915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(1.6289564) q[2];
sx q[2];
rz(-1.460961) q[2];
sx q[2];
rz(0.030208781) q[2];
rz(-0.31717832) q[3];
sx q[3];
rz(-0.89430292) q[3];
sx q[3];
rz(3.111307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
