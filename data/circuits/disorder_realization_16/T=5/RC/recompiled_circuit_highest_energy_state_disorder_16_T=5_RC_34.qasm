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
rz(-1.8104115) q[0];
sx q[0];
rz(-1.6348732) q[0];
sx q[0];
rz(2.9659086) q[0];
rz(-2.076258) q[1];
sx q[1];
rz(-1.330436) q[1];
sx q[1];
rz(0.65520823) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86533812) q[0];
sx q[0];
rz(-2.0575581) q[0];
sx q[0];
rz(-0.36051644) q[0];
rz(-pi) q[1];
rz(-1.6263481) q[2];
sx q[2];
rz(-1.4719988) q[2];
sx q[2];
rz(0.21888079) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53319327) q[1];
sx q[1];
rz(-3.0113868) q[1];
sx q[1];
rz(0.59334468) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1813222) q[3];
sx q[3];
rz(-1.8132079) q[3];
sx q[3];
rz(1.2038972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6931927) q[2];
sx q[2];
rz(-3.0312263) q[2];
sx q[2];
rz(2.8006862) q[2];
rz(0.77346268) q[3];
sx q[3];
rz(-2.1229459) q[3];
sx q[3];
rz(-3.036518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4648436) q[0];
sx q[0];
rz(-0.28443795) q[0];
sx q[0];
rz(-0.25962096) q[0];
rz(-2.3530841) q[1];
sx q[1];
rz(-2.2854243) q[1];
sx q[1];
rz(1.2451008) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21331025) q[0];
sx q[0];
rz(-1.6716372) q[0];
sx q[0];
rz(-2.1891602) q[0];
rz(1.4711667) q[2];
sx q[2];
rz(-1.2567695) q[2];
sx q[2];
rz(0.043836029) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3281365) q[1];
sx q[1];
rz(-2.3225647) q[1];
sx q[1];
rz(-2.0585895) q[1];
rz(-2.8471242) q[3];
sx q[3];
rz(-0.52669062) q[3];
sx q[3];
rz(1.9484438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.239324) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(-2.2321205) q[2];
rz(-2.4746312) q[3];
sx q[3];
rz(-1.704155) q[3];
sx q[3];
rz(0.10233574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7143836) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(-1.3877731) q[0];
rz(0.99675238) q[1];
sx q[1];
rz(-0.69147384) q[1];
sx q[1];
rz(0.48666418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62138689) q[0];
sx q[0];
rz(-2.7059816) q[0];
sx q[0];
rz(1.9389047) q[0];
rz(-pi) q[1];
rz(-0.66101199) q[2];
sx q[2];
rz(-2.1901988) q[2];
sx q[2];
rz(0.41601478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6909161) q[1];
sx q[1];
rz(-1.449739) q[1];
sx q[1];
rz(2.6208505) q[1];
rz(-pi) q[2];
rz(2.9938042) q[3];
sx q[3];
rz(-0.94919616) q[3];
sx q[3];
rz(-3.0382267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3643058) q[2];
sx q[2];
rz(-1.2591668) q[2];
sx q[2];
rz(0.49059179) q[2];
rz(-2.3250438) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(1.0402927) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33646026) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(2.1601987) q[0];
rz(-0.40311748) q[1];
sx q[1];
rz(-0.50794452) q[1];
sx q[1];
rz(-0.53781646) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92141672) q[0];
sx q[0];
rz(-2.8517637) q[0];
sx q[0];
rz(-2.1191862) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0079285) q[2];
sx q[2];
rz(-0.85459177) q[2];
sx q[2];
rz(1.4203526) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0723426) q[1];
sx q[1];
rz(-1.9543271) q[1];
sx q[1];
rz(1.0479397) q[1];
x q[2];
rz(1.1868401) q[3];
sx q[3];
rz(-2.4297485) q[3];
sx q[3];
rz(1.8432338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2526523) q[2];
sx q[2];
rz(-0.50709587) q[2];
sx q[2];
rz(1.0520774) q[2];
rz(-0.55225152) q[3];
sx q[3];
rz(-1.4058607) q[3];
sx q[3];
rz(-1.2205869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030516142) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(3.1040177) q[0];
rz(0.75282085) q[1];
sx q[1];
rz(-1.564874) q[1];
sx q[1];
rz(-0.084223025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1128453) q[0];
sx q[0];
rz(-1.0431093) q[0];
sx q[0];
rz(2.8535088) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2758314) q[2];
sx q[2];
rz(-1.2216481) q[2];
sx q[2];
rz(-0.63284208) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3854691) q[1];
sx q[1];
rz(-1.1471355) q[1];
sx q[1];
rz(-3.0961354) q[1];
rz(-pi) q[2];
rz(-1.8495848) q[3];
sx q[3];
rz(-2.9745347) q[3];
sx q[3];
rz(0.63705772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61421824) q[2];
sx q[2];
rz(-0.97443333) q[2];
sx q[2];
rz(-2.7091889) q[2];
rz(-1.4934941) q[3];
sx q[3];
rz(-1.3727539) q[3];
sx q[3];
rz(2.4026332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768141) q[0];
sx q[0];
rz(-1.5274763) q[0];
sx q[0];
rz(2.3798808) q[0];
rz(-1.0234045) q[1];
sx q[1];
rz(-0.49097148) q[1];
sx q[1];
rz(0.34058079) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0878076) q[0];
sx q[0];
rz(-1.4820966) q[0];
sx q[0];
rz(2.3596463) q[0];
rz(-2.5804065) q[2];
sx q[2];
rz(-1.7034987) q[2];
sx q[2];
rz(2.7787839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.611515) q[1];
sx q[1];
rz(-2.7002091) q[1];
sx q[1];
rz(-2.2468021) q[1];
x q[2];
rz(-0.16400985) q[3];
sx q[3];
rz(-0.7571547) q[3];
sx q[3];
rz(-1.8036158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1731825) q[2];
sx q[2];
rz(-2.0081655) q[2];
sx q[2];
rz(-0.75562149) q[2];
rz(0.38718811) q[3];
sx q[3];
rz(-1.3914934) q[3];
sx q[3];
rz(0.79425991) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8887535) q[0];
sx q[0];
rz(-3.0536953) q[0];
sx q[0];
rz(-1.2095691) q[0];
rz(-1.5225438) q[1];
sx q[1];
rz(-2.3698898) q[1];
sx q[1];
rz(-1.3873772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6722058) q[0];
sx q[0];
rz(-1.2325365) q[0];
sx q[0];
rz(-2.8126008) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31529324) q[2];
sx q[2];
rz(-1.0206501) q[2];
sx q[2];
rz(1.2619051) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55303326) q[1];
sx q[1];
rz(-0.84399763) q[1];
sx q[1];
rz(0.56808205) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.558616) q[3];
sx q[3];
rz(-2.600353) q[3];
sx q[3];
rz(-1.4547294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1986971) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(-0.4100619) q[2];
rz(-0.53747082) q[3];
sx q[3];
rz(-1.042807) q[3];
sx q[3];
rz(1.0679831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10802565) q[0];
sx q[0];
rz(-1.9060059) q[0];
sx q[0];
rz(1.0066475) q[0];
rz(2.0058696) q[1];
sx q[1];
rz(-0.4907116) q[1];
sx q[1];
rz(-0.27166414) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0245117) q[0];
sx q[0];
rz(-2.4371464) q[0];
sx q[0];
rz(-0.68450494) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92046787) q[2];
sx q[2];
rz(-2.1365385) q[2];
sx q[2];
rz(2.5793902) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9342088) q[1];
sx q[1];
rz(-2.1308239) q[1];
sx q[1];
rz(-0.70785849) q[1];
rz(-pi) q[2];
rz(2.397698) q[3];
sx q[3];
rz(-1.7366323) q[3];
sx q[3];
rz(1.3506571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2037619) q[2];
sx q[2];
rz(-2.2088642) q[2];
sx q[2];
rz(2.0358098) q[2];
rz(-1.3984937) q[3];
sx q[3];
rz(-2.1574557) q[3];
sx q[3];
rz(-2.2891146) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6451013) q[0];
sx q[0];
rz(-0.95106769) q[0];
sx q[0];
rz(-0.35354653) q[0];
rz(1.722466) q[1];
sx q[1];
rz(-1.6058233) q[1];
sx q[1];
rz(-3.0790192) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69819151) q[0];
sx q[0];
rz(-1.0834604) q[0];
sx q[0];
rz(0.62701012) q[0];
x q[1];
rz(2.563971) q[2];
sx q[2];
rz(-0.84813877) q[2];
sx q[2];
rz(1.3179145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8282916) q[1];
sx q[1];
rz(-2.5646659) q[1];
sx q[1];
rz(-2.5531656) q[1];
rz(-pi) q[2];
rz(-1.6639641) q[3];
sx q[3];
rz(-0.9103295) q[3];
sx q[3];
rz(2.4371393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.320437) q[2];
sx q[2];
rz(-1.7675567) q[2];
sx q[2];
rz(2.721411) q[2];
rz(-2.7968416) q[3];
sx q[3];
rz(-2.3638066) q[3];
sx q[3];
rz(0.93728089) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47278136) q[0];
sx q[0];
rz(-2.9720699) q[0];
sx q[0];
rz(1.8602759) q[0];
rz(1.5877113) q[1];
sx q[1];
rz(-2.3322767) q[1];
sx q[1];
rz(-0.70714998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4627461) q[0];
sx q[0];
rz(-0.94043676) q[0];
sx q[0];
rz(3.1267089) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98773445) q[2];
sx q[2];
rz(-1.1942004) q[2];
sx q[2];
rz(1.131112) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8445804) q[1];
sx q[1];
rz(-1.2606916) q[1];
sx q[1];
rz(-1.3363786) q[1];
rz(-pi) q[2];
rz(2.5159734) q[3];
sx q[3];
rz(-1.0418001) q[3];
sx q[3];
rz(2.4629018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3215434) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(0.496544) q[2];
rz(-1.2921565) q[3];
sx q[3];
rz(-0.29785952) q[3];
sx q[3];
rz(-2.7636102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9847066) q[0];
sx q[0];
rz(-2.8732185) q[0];
sx q[0];
rz(2.1068841) q[0];
rz(-2.5490419) q[1];
sx q[1];
rz(-2.4609346) q[1];
sx q[1];
rz(1.5998283) q[1];
rz(-1.8027923) q[2];
sx q[2];
rz(-1.6101888) q[2];
sx q[2];
rz(-2.0801795) q[2];
rz(1.3453398) q[3];
sx q[3];
rz(-2.1626948) q[3];
sx q[3];
rz(1.8656262) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
