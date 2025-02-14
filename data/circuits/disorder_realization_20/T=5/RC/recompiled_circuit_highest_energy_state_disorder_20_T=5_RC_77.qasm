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
rz(-2.9450077) q[0];
sx q[0];
rz(-0.44214806) q[0];
sx q[0];
rz(-2.2801939) q[0];
rz(1.5098894) q[1];
sx q[1];
rz(1.3246526) q[1];
sx q[1];
rz(9.4603705) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7073878) q[0];
sx q[0];
rz(-1.905283) q[0];
sx q[0];
rz(1.7777966) q[0];
x q[1];
rz(-0.33195095) q[2];
sx q[2];
rz(-0.32302007) q[2];
sx q[2];
rz(-0.92701757) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1110927) q[1];
sx q[1];
rz(-0.72734088) q[1];
sx q[1];
rz(0.05383454) q[1];
rz(-pi) q[2];
rz(2.2507812) q[3];
sx q[3];
rz(-1.5568131) q[3];
sx q[3];
rz(2.3774862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.55277905) q[2];
sx q[2];
rz(-1.8401044) q[2];
sx q[2];
rz(-2.0471052) q[2];
rz(-1.5158481) q[3];
sx q[3];
rz(-0.87259126) q[3];
sx q[3];
rz(-2.998013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0933519) q[0];
sx q[0];
rz(-0.99026647) q[0];
sx q[0];
rz(-2.7675203) q[0];
rz(-2.2834942) q[1];
sx q[1];
rz(-2.2007807) q[1];
sx q[1];
rz(-1.0757974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0556465) q[0];
sx q[0];
rz(-1.4688562) q[0];
sx q[0];
rz(1.6744958) q[0];
rz(-pi) q[1];
rz(-2.5941761) q[2];
sx q[2];
rz(-1.7760065) q[2];
sx q[2];
rz(2.2785435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75446426) q[1];
sx q[1];
rz(-2.7242047) q[1];
sx q[1];
rz(2.8512958) q[1];
rz(-0.76794736) q[3];
sx q[3];
rz(-2.30184) q[3];
sx q[3];
rz(-2.6382382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51521236) q[2];
sx q[2];
rz(-2.7084646) q[2];
sx q[2];
rz(-1.0832146) q[2];
rz(1.8488688) q[3];
sx q[3];
rz(-1.1336528) q[3];
sx q[3];
rz(2.2129272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2556297) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(-0.5249002) q[0];
rz(1.6820172) q[1];
sx q[1];
rz(-2.4039905) q[1];
sx q[1];
rz(2.9958013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6896714) q[0];
sx q[0];
rz(-1.4654241) q[0];
sx q[0];
rz(1.4793878) q[0];
x q[1];
rz(-0.96876933) q[2];
sx q[2];
rz(-2.0240236) q[2];
sx q[2];
rz(-0.22039686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.98602146) q[1];
sx q[1];
rz(-2.0029481) q[1];
sx q[1];
rz(-1.6701103) q[1];
rz(-pi) q[2];
rz(-0.94696857) q[3];
sx q[3];
rz(-0.83443975) q[3];
sx q[3];
rz(-0.26716993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8968481) q[2];
sx q[2];
rz(-1.1620099) q[2];
sx q[2];
rz(-2.9761918) q[2];
rz(0.8616972) q[3];
sx q[3];
rz(-2.4498037) q[3];
sx q[3];
rz(0.19680463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8621314) q[0];
sx q[0];
rz(-2.617351) q[0];
sx q[0];
rz(-0.72823802) q[0];
rz(2.8184452) q[1];
sx q[1];
rz(-1.0341945) q[1];
sx q[1];
rz(-1.039215) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208505) q[0];
sx q[0];
rz(-0.59950638) q[0];
sx q[0];
rz(-2.3248424) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2065998) q[2];
sx q[2];
rz(-1.3424917) q[2];
sx q[2];
rz(-1.993597) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8700712) q[1];
sx q[1];
rz(-1.0876552) q[1];
sx q[1];
rz(1.4663694) q[1];
rz(0.69833243) q[3];
sx q[3];
rz(-0.94925967) q[3];
sx q[3];
rz(-0.96635287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6497961) q[2];
sx q[2];
rz(-2.0911262) q[2];
sx q[2];
rz(0.51898471) q[2];
rz(-2.0670048) q[3];
sx q[3];
rz(-1.2858177) q[3];
sx q[3];
rz(-2.6463267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.062926) q[0];
sx q[0];
rz(-0.92200297) q[0];
sx q[0];
rz(2.9700188) q[0];
rz(1.0942787) q[1];
sx q[1];
rz(-1.3010052) q[1];
sx q[1];
rz(0.23460728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6162524) q[0];
sx q[0];
rz(-1.293129) q[0];
sx q[0];
rz(-3.0086631) q[0];
rz(2.634356) q[2];
sx q[2];
rz(-1.5248858) q[2];
sx q[2];
rz(-3.1243589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.34688309) q[1];
sx q[1];
rz(-0.88635072) q[1];
sx q[1];
rz(2.3331932) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9453641) q[3];
sx q[3];
rz(-0.8959594) q[3];
sx q[3];
rz(-0.57317153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91308633) q[2];
sx q[2];
rz(-1.9501016) q[2];
sx q[2];
rz(1.4027493) q[2];
rz(2.5168822) q[3];
sx q[3];
rz(-0.96840817) q[3];
sx q[3];
rz(1.7984084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45873555) q[0];
sx q[0];
rz(-0.0077489297) q[0];
sx q[0];
rz(-2.8140581) q[0];
rz(0.028118357) q[1];
sx q[1];
rz(-1.997812) q[1];
sx q[1];
rz(-2.1498674) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018143749) q[0];
sx q[0];
rz(-0.94052343) q[0];
sx q[0];
rz(2.2176377) q[0];
rz(-0.2260416) q[2];
sx q[2];
rz(-1.4648572) q[2];
sx q[2];
rz(-0.49913479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2497445) q[1];
sx q[1];
rz(-1.6127024) q[1];
sx q[1];
rz(-0.81467198) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47073165) q[3];
sx q[3];
rz(-1.5464029) q[3];
sx q[3];
rz(-0.40917802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.907054) q[2];
sx q[2];
rz(-1.3776642) q[2];
sx q[2];
rz(-1.6592337) q[2];
rz(2.0756857) q[3];
sx q[3];
rz(-1.1803455) q[3];
sx q[3];
rz(-2.0445686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3858353) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(0.29543153) q[0];
rz(2.9040728) q[1];
sx q[1];
rz(-0.93370456) q[1];
sx q[1];
rz(0.74737731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45493653) q[0];
sx q[0];
rz(-2.1224408) q[0];
sx q[0];
rz(0.96740361) q[0];
rz(-2.5146444) q[2];
sx q[2];
rz(-2.0184506) q[2];
sx q[2];
rz(2.3989776) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.69773) q[1];
sx q[1];
rz(-1.2011114) q[1];
sx q[1];
rz(-2.1013942) q[1];
rz(-pi) q[2];
rz(0.5434955) q[3];
sx q[3];
rz(-1.8955909) q[3];
sx q[3];
rz(2.3372758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4733009) q[2];
sx q[2];
rz(-1.1857727) q[2];
sx q[2];
rz(0.2549003) q[2];
rz(0.21236803) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(-3.1406241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2905228) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(3.0176924) q[0];
rz(-0.87521416) q[1];
sx q[1];
rz(-0.83299914) q[1];
sx q[1];
rz(-2.5828054) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0103115) q[0];
sx q[0];
rz(-2.0925267) q[0];
sx q[0];
rz(-0.2894131) q[0];
x q[1];
rz(2.8577639) q[2];
sx q[2];
rz(-2.0528194) q[2];
sx q[2];
rz(1.986077) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9732157) q[1];
sx q[1];
rz(-1.3011908) q[1];
sx q[1];
rz(2.5942179) q[1];
rz(-2.111192) q[3];
sx q[3];
rz(-2.4321788) q[3];
sx q[3];
rz(-1.7075001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.045804068) q[2];
sx q[2];
rz(-2.6498821) q[2];
sx q[2];
rz(2.4208505) q[2];
rz(2.7785684) q[3];
sx q[3];
rz(-1.9178773) q[3];
sx q[3];
rz(1.5117517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1849798) q[0];
sx q[0];
rz(-2.9473801) q[0];
sx q[0];
rz(-0.089056253) q[0];
rz(2.7748499) q[1];
sx q[1];
rz(-2.2373503) q[1];
sx q[1];
rz(-1.6820224) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3066481) q[0];
sx q[0];
rz(-2.1900926) q[0];
sx q[0];
rz(0.28138103) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4602997) q[2];
sx q[2];
rz(-0.80901399) q[2];
sx q[2];
rz(-1.9495277) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0017346026) q[1];
sx q[1];
rz(-2.325104) q[1];
sx q[1];
rz(1.0553318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60114558) q[3];
sx q[3];
rz(-2.450305) q[3];
sx q[3];
rz(2.9959681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9350932) q[2];
sx q[2];
rz(-1.3672028) q[2];
sx q[2];
rz(1.8915141) q[2];
rz(-1.2262723) q[3];
sx q[3];
rz(-0.63392249) q[3];
sx q[3];
rz(0.63898501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6298237) q[0];
sx q[0];
rz(-2.0084232) q[0];
sx q[0];
rz(1.1400219) q[0];
rz(2.5584768) q[1];
sx q[1];
rz(-1.6551599) q[1];
sx q[1];
rz(0.66782943) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44396469) q[0];
sx q[0];
rz(-2.2831342) q[0];
sx q[0];
rz(1.1348073) q[0];
x q[1];
rz(0.23867757) q[2];
sx q[2];
rz(-0.88310034) q[2];
sx q[2];
rz(-0.65438813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9676303) q[1];
sx q[1];
rz(-0.47904992) q[1];
sx q[1];
rz(-1.1381989) q[1];
rz(-0.79270122) q[3];
sx q[3];
rz(-1.2374733) q[3];
sx q[3];
rz(1.9938716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84393152) q[2];
sx q[2];
rz(-2.6466978) q[2];
sx q[2];
rz(0.75221357) q[2];
rz(-3.0609868) q[3];
sx q[3];
rz(-1.7919431) q[3];
sx q[3];
rz(-0.54752553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053631) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(2.2930131) q[1];
sx q[1];
rz(-2.4129557) q[1];
sx q[1];
rz(1.6853263) q[1];
rz(2.3586629) q[2];
sx q[2];
rz(-2.30884) q[2];
sx q[2];
rz(-2.5306551) q[2];
rz(2.4918588) q[3];
sx q[3];
rz(-0.55745535) q[3];
sx q[3];
rz(-1.0563323) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
