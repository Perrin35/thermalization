OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2562113) q[0];
sx q[0];
rz(5.5698759) q[0];
sx q[0];
rz(8.2110693) q[0];
rz(1.775939) q[1];
sx q[1];
rz(-0.27888271) q[1];
sx q[1];
rz(-2.9409148) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0943692) q[0];
sx q[0];
rz(-0.46877623) q[0];
sx q[0];
rz(-1.3865406) q[0];
rz(-0.63849475) q[2];
sx q[2];
rz(-2.8732277) q[2];
sx q[2];
rz(-0.61682781) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38814784) q[1];
sx q[1];
rz(-2.8737443) q[1];
sx q[1];
rz(-0.74117383) q[1];
x q[2];
rz(-1.3986271) q[3];
sx q[3];
rz(-1.2667613) q[3];
sx q[3];
rz(-1.1393169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.063244907) q[2];
sx q[2];
rz(-1.1996317) q[2];
sx q[2];
rz(-2.4285748) q[2];
rz(1.2837422) q[3];
sx q[3];
rz(-2.4132437) q[3];
sx q[3];
rz(-0.89948765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17083183) q[0];
sx q[0];
rz(-1.4749227) q[0];
sx q[0];
rz(-1.4738039) q[0];
rz(-0.57227349) q[1];
sx q[1];
rz(-1.0621366) q[1];
sx q[1];
rz(-1.3776113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3085605) q[0];
sx q[0];
rz(-1.8695117) q[0];
sx q[0];
rz(2.0479982) q[0];
rz(2.6102561) q[2];
sx q[2];
rz(-1.9730933) q[2];
sx q[2];
rz(-0.50237331) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2050537) q[1];
sx q[1];
rz(-1.8636101) q[1];
sx q[1];
rz(-1.9163301) q[1];
rz(0.99239852) q[3];
sx q[3];
rz(-2.7790894) q[3];
sx q[3];
rz(0.9622919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1231692) q[2];
sx q[2];
rz(-1.5095242) q[2];
sx q[2];
rz(1.8040166) q[2];
rz(0.34559524) q[3];
sx q[3];
rz(-2.5497422) q[3];
sx q[3];
rz(0.93446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086394101) q[0];
sx q[0];
rz(-0.18562695) q[0];
sx q[0];
rz(-0.58919543) q[0];
rz(-0.20801726) q[1];
sx q[1];
rz(-0.55744019) q[1];
sx q[1];
rz(0.20588188) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7939121) q[0];
sx q[0];
rz(-0.97647515) q[0];
sx q[0];
rz(-1.1403313) q[0];
rz(-2.6462502) q[2];
sx q[2];
rz(-2.2215507) q[2];
sx q[2];
rz(-2.9913354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66135397) q[1];
sx q[1];
rz(-1.2671372) q[1];
sx q[1];
rz(-2.420029) q[1];
rz(-0.08783665) q[3];
sx q[3];
rz(-1.1655775) q[3];
sx q[3];
rz(2.3072568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0899293) q[2];
sx q[2];
rz(-2.0970924) q[2];
sx q[2];
rz(0.54452407) q[2];
rz(0.45474592) q[3];
sx q[3];
rz(-2.5179722) q[3];
sx q[3];
rz(0.0039984306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74314463) q[0];
sx q[0];
rz(-2.5938617) q[0];
sx q[0];
rz(-0.69993436) q[0];
rz(-1.3088538) q[1];
sx q[1];
rz(-2.0715641) q[1];
sx q[1];
rz(-1.4189789) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5154079) q[0];
sx q[0];
rz(-1.1061449) q[0];
sx q[0];
rz(0.22494577) q[0];
x q[1];
rz(1.5350902) q[2];
sx q[2];
rz(-0.49884819) q[2];
sx q[2];
rz(0.6998261) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0154889) q[1];
sx q[1];
rz(-2.2601193) q[1];
sx q[1];
rz(1.2587121) q[1];
rz(-pi) q[2];
rz(-0.14601018) q[3];
sx q[3];
rz(-1.8838804) q[3];
sx q[3];
rz(-2.2816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2122638) q[2];
sx q[2];
rz(-0.45390359) q[2];
sx q[2];
rz(-1.7960499) q[2];
rz(-2.4228607) q[3];
sx q[3];
rz(-0.74145442) q[3];
sx q[3];
rz(-2.45347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5162002) q[0];
sx q[0];
rz(-1.1855519) q[0];
sx q[0];
rz(0.54310435) q[0];
rz(2.886046) q[1];
sx q[1];
rz(-2.6592022) q[1];
sx q[1];
rz(1.3822752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4323174) q[0];
sx q[0];
rz(-2.0943644) q[0];
sx q[0];
rz(-1.104959) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7167983) q[2];
sx q[2];
rz(-0.80933648) q[2];
sx q[2];
rz(-2.3392364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4529323) q[1];
sx q[1];
rz(-2.1218461) q[1];
sx q[1];
rz(0.30639415) q[1];
rz(-pi) q[2];
rz(2.7409389) q[3];
sx q[3];
rz(-1.7285498) q[3];
sx q[3];
rz(0.98689134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1725267) q[2];
sx q[2];
rz(-1.669599) q[2];
sx q[2];
rz(-2.772061) q[2];
rz(-1.8116123) q[3];
sx q[3];
rz(-0.79052916) q[3];
sx q[3];
rz(1.7867521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0238277) q[0];
sx q[0];
rz(-2.8067639) q[0];
sx q[0];
rz(-0.086300015) q[0];
rz(-1.49508) q[1];
sx q[1];
rz(-1.1777271) q[1];
sx q[1];
rz(0.42974791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3783898) q[0];
sx q[0];
rz(-2.9400005) q[0];
sx q[0];
rz(1.2132573) q[0];
rz(-pi) q[1];
rz(0.65799539) q[2];
sx q[2];
rz(-1.3086196) q[2];
sx q[2];
rz(-3.0397751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6657998) q[1];
sx q[1];
rz(-2.0415039) q[1];
sx q[1];
rz(-1.0757955) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67062258) q[3];
sx q[3];
rz(-1.792893) q[3];
sx q[3];
rz(2.5012453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0325844) q[2];
sx q[2];
rz(-2.0592368) q[2];
sx q[2];
rz(-1.2004131) q[2];
rz(0.16610185) q[3];
sx q[3];
rz(-0.76986543) q[3];
sx q[3];
rz(-2.2782245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2849041) q[0];
sx q[0];
rz(-0.76304522) q[0];
sx q[0];
rz(2.3175008) q[0];
rz(1.9169982) q[1];
sx q[1];
rz(-1.9173744) q[1];
sx q[1];
rz(1.0138938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4606095) q[0];
sx q[0];
rz(-1.8797726) q[0];
sx q[0];
rz(0.81134422) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5297671) q[2];
sx q[2];
rz(-2.6227747) q[2];
sx q[2];
rz(-0.045836115) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2353603) q[1];
sx q[1];
rz(-2.4543833) q[1];
sx q[1];
rz(-0.63628025) q[1];
rz(-pi) q[2];
rz(-0.9634094) q[3];
sx q[3];
rz(-0.62650354) q[3];
sx q[3];
rz(2.1914848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27807221) q[2];
sx q[2];
rz(-2.4188228) q[2];
sx q[2];
rz(-1.2449167) q[2];
rz(-2.909929) q[3];
sx q[3];
rz(-1.04117) q[3];
sx q[3];
rz(0.1483354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0949377) q[0];
sx q[0];
rz(-1.2060839) q[0];
sx q[0];
rz(0.99223247) q[0];
rz(-0.16920432) q[1];
sx q[1];
rz(-0.84183401) q[1];
sx q[1];
rz(-1.4134891) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0071825562) q[0];
sx q[0];
rz(-2.5408486) q[0];
sx q[0];
rz(-0.13237615) q[0];
rz(-1.0186853) q[2];
sx q[2];
rz(-1.8183961) q[2];
sx q[2];
rz(-1.1820861) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0553356) q[1];
sx q[1];
rz(-0.29976832) q[1];
sx q[1];
rz(-0.10792984) q[1];
rz(-pi) q[2];
rz(2.9417073) q[3];
sx q[3];
rz(-1.098714) q[3];
sx q[3];
rz(2.1560644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32470545) q[2];
sx q[2];
rz(-1.4810666) q[2];
sx q[2];
rz(-1.7186349) q[2];
rz(0.11296806) q[3];
sx q[3];
rz(-0.26765099) q[3];
sx q[3];
rz(0.62937361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896352) q[0];
sx q[0];
rz(-1.4653787) q[0];
sx q[0];
rz(-2.5919609) q[0];
rz(-2.5426087) q[1];
sx q[1];
rz(-2.2689029) q[1];
sx q[1];
rz(-1.6302861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0881168) q[0];
sx q[0];
rz(-1.5673182) q[0];
sx q[0];
rz(-0.17668488) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51038607) q[2];
sx q[2];
rz(-0.56737075) q[2];
sx q[2];
rz(0.20287831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28569878) q[1];
sx q[1];
rz(-2.0715069) q[1];
sx q[1];
rz(0.46127528) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6069943) q[3];
sx q[3];
rz(-0.55639297) q[3];
sx q[3];
rz(-2.8884187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0281684) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(-2.8070731) q[2];
rz(0.64535514) q[3];
sx q[3];
rz(-2.1367475) q[3];
sx q[3];
rz(0.2373124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065011218) q[0];
sx q[0];
rz(-0.33877057) q[0];
sx q[0];
rz(0.6231128) q[0];
rz(0.30162853) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(-1.041144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17693263) q[0];
sx q[0];
rz(-0.76632351) q[0];
sx q[0];
rz(-2.4936952) q[0];
x q[1];
rz(0.1847965) q[2];
sx q[2];
rz(-0.46773887) q[2];
sx q[2];
rz(-2.6370905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1438801) q[1];
sx q[1];
rz(-1.5399944) q[1];
sx q[1];
rz(0.7672337) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0431533) q[3];
sx q[3];
rz(-2.7343547) q[3];
sx q[3];
rz(-1.094363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4447896) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(0.94318548) q[2];
rz(1.4140363) q[3];
sx q[3];
rz(-0.90641886) q[3];
sx q[3];
rz(0.79677719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74054756) q[0];
sx q[0];
rz(-2.7475806) q[0];
sx q[0];
rz(-0.077234118) q[0];
rz(-0.41863353) q[1];
sx q[1];
rz(-1.5593465) q[1];
sx q[1];
rz(-2.9895463) q[1];
rz(-0.05337333) q[2];
sx q[2];
rz(-2.5659701) q[2];
sx q[2];
rz(0.26502668) q[2];
rz(1.2996243) q[3];
sx q[3];
rz(-1.9883755) q[3];
sx q[3];
rz(-1.1379013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
