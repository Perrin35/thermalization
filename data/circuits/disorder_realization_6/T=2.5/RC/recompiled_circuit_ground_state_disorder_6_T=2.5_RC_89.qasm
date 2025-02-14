OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4226469) q[0];
sx q[0];
rz(-0.5889686) q[0];
sx q[0];
rz(0.11864057) q[0];
rz(-0.9912107) q[1];
sx q[1];
rz(4.1832357) q[1];
sx q[1];
rz(15.19419) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0002961) q[0];
sx q[0];
rz(-2.0679255) q[0];
sx q[0];
rz(0.5750467) q[0];
rz(-pi) q[1];
rz(1.5085601) q[2];
sx q[2];
rz(-1.0966435) q[2];
sx q[2];
rz(1.7218423) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.895925) q[1];
sx q[1];
rz(-1.8918719) q[1];
sx q[1];
rz(1.7968057) q[1];
x q[2];
rz(-0.59207497) q[3];
sx q[3];
rz(-1.5993803) q[3];
sx q[3];
rz(-2.0297272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4673246) q[2];
sx q[2];
rz(-1.5200204) q[2];
sx q[2];
rz(-0.33217126) q[2];
rz(0.25003555) q[3];
sx q[3];
rz(-1.2269521) q[3];
sx q[3];
rz(1.2927607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2501204) q[0];
sx q[0];
rz(-1.8880867) q[0];
sx q[0];
rz(2.6943595) q[0];
rz(1.0294634) q[1];
sx q[1];
rz(-0.58164683) q[1];
sx q[1];
rz(0.10118016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44199884) q[0];
sx q[0];
rz(-1.170982) q[0];
sx q[0];
rz(0.48161048) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1498141) q[2];
sx q[2];
rz(-1.4027032) q[2];
sx q[2];
rz(1.6630295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.356368) q[1];
sx q[1];
rz(-1.984298) q[1];
sx q[1];
rz(1.2797536) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4199791) q[3];
sx q[3];
rz(-2.4972417) q[3];
sx q[3];
rz(2.2805813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0589361) q[2];
sx q[2];
rz(-0.75581789) q[2];
sx q[2];
rz(1.5400881) q[2];
rz(-0.61795175) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(1.1624973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5471632) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(-2.6255703) q[0];
rz(1.0524606) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(1.9693536) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7965423) q[0];
sx q[0];
rz(-1.5818514) q[0];
sx q[0];
rz(0.42406908) q[0];
x q[1];
rz(-0.86031057) q[2];
sx q[2];
rz(-1.6261618) q[2];
sx q[2];
rz(0.64337984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8015299) q[1];
sx q[1];
rz(-1.0822359) q[1];
sx q[1];
rz(-2.9304781) q[1];
x q[2];
rz(-3.1276305) q[3];
sx q[3];
rz(-2.3413071) q[3];
sx q[3];
rz(-1.6228907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8664794) q[2];
sx q[2];
rz(-1.2629513) q[2];
sx q[2];
rz(0.67683721) q[2];
rz(2.6721241) q[3];
sx q[3];
rz(-2.2468086) q[3];
sx q[3];
rz(1.9278056) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6770099) q[0];
sx q[0];
rz(-2.1362169) q[0];
sx q[0];
rz(-1.1997892) q[0];
rz(2.5489573) q[1];
sx q[1];
rz(-1.0721782) q[1];
sx q[1];
rz(-1.6823654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98226794) q[0];
sx q[0];
rz(-1.3749116) q[0];
sx q[0];
rz(-2.5513493) q[0];
x q[1];
rz(0.24942579) q[2];
sx q[2];
rz(-0.81455961) q[2];
sx q[2];
rz(-0.047175353) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4479936) q[1];
sx q[1];
rz(-1.3136615) q[1];
sx q[1];
rz(0.93235459) q[1];
rz(0.45444416) q[3];
sx q[3];
rz(-2.4049149) q[3];
sx q[3];
rz(-3.122913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49372855) q[2];
sx q[2];
rz(-1.0961327) q[2];
sx q[2];
rz(0.68191051) q[2];
rz(-2.72825) q[3];
sx q[3];
rz(-1.5324493) q[3];
sx q[3];
rz(1.9857064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2979564) q[0];
sx q[0];
rz(-1.4763259) q[0];
sx q[0];
rz(-2.8635039) q[0];
rz(-0.9043215) q[1];
sx q[1];
rz(-1.6878637) q[1];
sx q[1];
rz(2.1542737) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7287059) q[0];
sx q[0];
rz(-0.97620076) q[0];
sx q[0];
rz(0.28810687) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1369517) q[2];
sx q[2];
rz(-1.5986575) q[2];
sx q[2];
rz(0.34649039) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4523068) q[1];
sx q[1];
rz(-2.747162) q[1];
sx q[1];
rz(1.4107262) q[1];
rz(-2.0102536) q[3];
sx q[3];
rz(-0.58142291) q[3];
sx q[3];
rz(1.9245337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.014293369) q[2];
sx q[2];
rz(-2.1743446) q[2];
sx q[2];
rz(-0.75931749) q[2];
rz(1.4826108) q[3];
sx q[3];
rz(-1.1967775) q[3];
sx q[3];
rz(-0.35879859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21869126) q[0];
sx q[0];
rz(-0.8194812) q[0];
sx q[0];
rz(-2.4500093) q[0];
rz(2.2870731) q[1];
sx q[1];
rz(-2.4311192) q[1];
sx q[1];
rz(-0.42133322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92625657) q[0];
sx q[0];
rz(-0.43784446) q[0];
sx q[0];
rz(-1.7329751) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3106903) q[2];
sx q[2];
rz(-2.75616) q[2];
sx q[2];
rz(-1.3192434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5798963) q[1];
sx q[1];
rz(-1.6547799) q[1];
sx q[1];
rz(-2.7799699) q[1];
x q[2];
rz(1.8438897) q[3];
sx q[3];
rz(-2.1583302) q[3];
sx q[3];
rz(-1.7949826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.93115807) q[2];
sx q[2];
rz(-1.4782108) q[2];
sx q[2];
rz(2.9388536) q[2];
rz(1.5984795) q[3];
sx q[3];
rz(-0.32035443) q[3];
sx q[3];
rz(1.6725484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8572674) q[0];
sx q[0];
rz(-1.7131282) q[0];
sx q[0];
rz(-0.39886928) q[0];
rz(1.9786037) q[1];
sx q[1];
rz(-2.3154924) q[1];
sx q[1];
rz(-0.2805447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60567509) q[0];
sx q[0];
rz(-1.1725764) q[0];
sx q[0];
rz(-1.1541973) q[0];
x q[1];
rz(1.0516099) q[2];
sx q[2];
rz(-2.3943442) q[2];
sx q[2];
rz(-2.0488536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7458492) q[1];
sx q[1];
rz(-1.1458995) q[1];
sx q[1];
rz(2.8060071) q[1];
x q[2];
rz(0.58275025) q[3];
sx q[3];
rz(-2.2556335) q[3];
sx q[3];
rz(-0.69151758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3205388) q[2];
sx q[2];
rz(-0.97964779) q[2];
sx q[2];
rz(1.4873571) q[2];
rz(2.1982543) q[3];
sx q[3];
rz(-2.2127559) q[3];
sx q[3];
rz(1.1613891) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6134375) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(-2.6392537) q[0];
rz(2.595064) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(-2.6796403) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97945178) q[0];
sx q[0];
rz(-0.58425036) q[0];
sx q[0];
rz(-2.1186078) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0937198) q[2];
sx q[2];
rz(-0.92000735) q[2];
sx q[2];
rz(1.6597468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89901224) q[1];
sx q[1];
rz(-1.1625152) q[1];
sx q[1];
rz(0.11870833) q[1];
x q[2];
rz(-0.53047385) q[3];
sx q[3];
rz(-1.5598205) q[3];
sx q[3];
rz(-0.61179841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.97303566) q[2];
sx q[2];
rz(-0.23739561) q[2];
sx q[2];
rz(0.17042223) q[2];
rz(-0.45525822) q[3];
sx q[3];
rz(-1.5452496) q[3];
sx q[3];
rz(-0.71793238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97487226) q[0];
sx q[0];
rz(-1.6396739) q[0];
sx q[0];
rz(-0.18044743) q[0];
rz(0.51668733) q[1];
sx q[1];
rz(-1.5770715) q[1];
sx q[1];
rz(0.43389854) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2045) q[0];
sx q[0];
rz(-2.4334416) q[0];
sx q[0];
rz(-2.2592179) q[0];
rz(-pi) q[1];
rz(2.216748) q[2];
sx q[2];
rz(-0.5415104) q[2];
sx q[2];
rz(-2.4157267) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5499159) q[1];
sx q[1];
rz(-2.2001451) q[1];
sx q[1];
rz(-2.9931328) q[1];
rz(-2.850344) q[3];
sx q[3];
rz(-2.2403803) q[3];
sx q[3];
rz(2.254829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58069289) q[2];
sx q[2];
rz(-0.43792024) q[2];
sx q[2];
rz(1.8221347) q[2];
rz(2.8699919) q[3];
sx q[3];
rz(-1.8235794) q[3];
sx q[3];
rz(-2.0297089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44529799) q[0];
sx q[0];
rz(-2.2309208) q[0];
sx q[0];
rz(1.4060422) q[0];
rz(-1.4259074) q[1];
sx q[1];
rz(-1.1049263) q[1];
sx q[1];
rz(1.2474733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96811622) q[0];
sx q[0];
rz(-2.0531468) q[0];
sx q[0];
rz(0.29598693) q[0];
x q[1];
rz(-0.054385238) q[2];
sx q[2];
rz(-1.453389) q[2];
sx q[2];
rz(-0.94737999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0669603) q[1];
sx q[1];
rz(-1.5299284) q[1];
sx q[1];
rz(2.7994179) q[1];
x q[2];
rz(-2.3214809) q[3];
sx q[3];
rz(-0.81104987) q[3];
sx q[3];
rz(0.91203989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71331104) q[2];
sx q[2];
rz(-0.2321299) q[2];
sx q[2];
rz(0.10078079) q[2];
rz(-0.0095602592) q[3];
sx q[3];
rz(-1.5848426) q[3];
sx q[3];
rz(-1.5338219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1495001) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(-0.59303444) q[1];
sx q[1];
rz(-1.4321764) q[1];
sx q[1];
rz(2.166116) q[1];
rz(-0.070746919) q[2];
sx q[2];
rz(-1.1732709) q[2];
sx q[2];
rz(2.7931961) q[2];
rz(-0.22105263) q[3];
sx q[3];
rz(-1.5536158) q[3];
sx q[3];
rz(2.0719846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
