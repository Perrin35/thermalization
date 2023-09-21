OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.68552652) q[0];
sx q[0];
rz(-2.752562) q[0];
sx q[0];
rz(0.88357893) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(-1.6844123) q[1];
sx q[1];
rz(-1.943346) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0966914) q[0];
sx q[0];
rz(-2.9400819) q[0];
sx q[0];
rz(0.50045307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3660907) q[2];
sx q[2];
rz(-1.8047793) q[2];
sx q[2];
rz(0.50228679) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7361684) q[1];
sx q[1];
rz(-0.6286469) q[1];
sx q[1];
rz(-0.58341649) q[1];
rz(-1.7415813) q[3];
sx q[3];
rz(-1.6163974) q[3];
sx q[3];
rz(-0.036269773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9464232) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(0.18134376) q[2];
rz(0.26120734) q[3];
sx q[3];
rz(-1.2657335) q[3];
sx q[3];
rz(-2.3852824) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(-0.38683495) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5418672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4138448) q[0];
sx q[0];
rz(-1.1607329) q[0];
sx q[0];
rz(-0.64322612) q[0];
rz(-1.8962757) q[2];
sx q[2];
rz(-2.6633334) q[2];
sx q[2];
rz(-1.9667728) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9215556) q[1];
sx q[1];
rz(-1.6065292) q[1];
sx q[1];
rz(-1.5608556) q[1];
rz(-pi) q[2];
rz(-2.6582791) q[3];
sx q[3];
rz(-1.1477074) q[3];
sx q[3];
rz(-1.2276633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5388422) q[2];
sx q[2];
rz(-0.87783146) q[2];
sx q[2];
rz(1.7269469) q[2];
rz(-0.85033068) q[3];
sx q[3];
rz(-2.7089705) q[3];
sx q[3];
rz(-1.458228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8787815) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(-0.61022726) q[0];
rz(1.2894851) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(-0.99197018) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2881644) q[0];
sx q[0];
rz(-0.58840226) q[0];
sx q[0];
rz(1.9073652) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8875071) q[2];
sx q[2];
rz(-2.457329) q[2];
sx q[2];
rz(2.4917045) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75533463) q[1];
sx q[1];
rz(-1.4812246) q[1];
sx q[1];
rz(2.7002525) q[1];
x q[2];
rz(-1.3258341) q[3];
sx q[3];
rz(-1.1253329) q[3];
sx q[3];
rz(-0.20416343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.38301864) q[2];
sx q[2];
rz(-0.16123161) q[2];
sx q[2];
rz(0.26829159) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-0.18850732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2878993) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(2.7365141) q[0];
rz(-0.69008094) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(-0.69782034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47397428) q[0];
sx q[0];
rz(-1.4809429) q[0];
sx q[0];
rz(-3.1180179) q[0];
rz(-pi) q[1];
rz(1.9525098) q[2];
sx q[2];
rz(-1.0587947) q[2];
sx q[2];
rz(-0.94101671) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3541672) q[1];
sx q[1];
rz(-2.0560871) q[1];
sx q[1];
rz(-0.20809681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9911733) q[3];
sx q[3];
rz(-1.2983409) q[3];
sx q[3];
rz(3.0416995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.5092124) q[2];
rz(-0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8835835) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(2.1160545) q[0];
rz(-2.569596) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(-0.62932032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6228018) q[0];
sx q[0];
rz(-0.51827058) q[0];
sx q[0];
rz(-0.63664125) q[0];
x q[1];
rz(2.2467381) q[2];
sx q[2];
rz(-1.7712777) q[2];
sx q[2];
rz(1.5467874) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3931261) q[1];
sx q[1];
rz(-1.6854291) q[1];
sx q[1];
rz(-0.51490358) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26341565) q[3];
sx q[3];
rz(-0.98589555) q[3];
sx q[3];
rz(-1.6057305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59262529) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(-0.34234753) q[2];
rz(-1.3458378) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489007) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(2.956399) q[0];
rz(-1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(-1.7746183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0104116) q[0];
sx q[0];
rz(-2.512012) q[0];
sx q[0];
rz(-2.7601943) q[0];
rz(1.9063437) q[2];
sx q[2];
rz(-2.3830072) q[2];
sx q[2];
rz(-0.89546613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1176599) q[1];
sx q[1];
rz(-1.4791094) q[1];
sx q[1];
rz(-1.4183284) q[1];
rz(1.1432511) q[3];
sx q[3];
rz(-2.098009) q[3];
sx q[3];
rz(-2.9851819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8967445) q[2];
sx q[2];
rz(-2.5881793) q[2];
sx q[2];
rz(0.883376) q[2];
rz(1.7287438) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(-0.81378716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8900523) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(-1.863377) q[0];
rz(3.1037519) q[1];
sx q[1];
rz(-0.81532878) q[1];
sx q[1];
rz(1.7657123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85526953) q[0];
sx q[0];
rz(-2.3648242) q[0];
sx q[0];
rz(-1.1458678) q[0];
x q[1];
rz(-1.4543578) q[2];
sx q[2];
rz(-2.6652626) q[2];
sx q[2];
rz(-1.6778698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.081825) q[1];
sx q[1];
rz(-1.4653112) q[1];
sx q[1];
rz(-3.0256773) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7405628) q[3];
sx q[3];
rz(-2.8068636) q[3];
sx q[3];
rz(1.7498121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4574796) q[2];
sx q[2];
rz(-2.423968) q[2];
sx q[2];
rz(-0.73105556) q[2];
rz(-3.030792) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59657997) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(-0.73356432) q[0];
rz(0.60797524) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(-0.2342934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83345862) q[0];
sx q[0];
rz(-1.3123543) q[0];
sx q[0];
rz(1.0779557) q[0];
rz(-pi) q[1];
rz(-2.2393164) q[2];
sx q[2];
rz(-2.1736645) q[2];
sx q[2];
rz(2.5831646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81208166) q[1];
sx q[1];
rz(-0.78204621) q[1];
sx q[1];
rz(-1.8409607) q[1];
rz(-pi) q[2];
rz(1.0844564) q[3];
sx q[3];
rz(-0.64389766) q[3];
sx q[3];
rz(-0.38734303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0155448) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(1.7162494) q[2];
rz(-1.6783293) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(-1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5995246) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(1.9482127) q[0];
rz(-1.880973) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(2.1967922) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017100447) q[0];
sx q[0];
rz(-1.955535) q[0];
sx q[0];
rz(2.4186633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1759042) q[2];
sx q[2];
rz(-1.8443622) q[2];
sx q[2];
rz(-0.045189518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24430844) q[1];
sx q[1];
rz(-2.8123887) q[1];
sx q[1];
rz(0.55397482) q[1];
rz(-pi) q[2];
rz(-1.2148083) q[3];
sx q[3];
rz(-1.755852) q[3];
sx q[3];
rz(-1.9389648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9251359) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(2.8533868) q[2];
rz(0.47973412) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(-2.2286041) q[0];
rz(0.37462014) q[1];
sx q[1];
rz(-1.7381784) q[1];
sx q[1];
rz(2.250681) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8091781) q[0];
sx q[0];
rz(-2.7190468) q[0];
sx q[0];
rz(-0.64230625) q[0];
rz(2.183379) q[2];
sx q[2];
rz(-1.805086) q[2];
sx q[2];
rz(2.2189552) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9526457) q[1];
sx q[1];
rz(-1.2332321) q[1];
sx q[1];
rz(1.360421) q[1];
x q[2];
rz(-2.9534146) q[3];
sx q[3];
rz(-1.6450226) q[3];
sx q[3];
rz(-1.0159514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23641071) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(0.50160828) q[2];
rz(1.2891399) q[3];
sx q[3];
rz(-1.6882378) q[3];
sx q[3];
rz(2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5065153) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(1.1322017) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(-2.6272527) q[2];
sx q[2];
rz(-2.8375576) q[2];
sx q[2];
rz(0.67630771) q[2];
rz(-1.6518456) q[3];
sx q[3];
rz(-0.54809082) q[3];
sx q[3];
rz(-2.5561668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
