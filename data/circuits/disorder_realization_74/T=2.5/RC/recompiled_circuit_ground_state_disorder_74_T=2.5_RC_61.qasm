OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.35535204) q[0];
sx q[0];
rz(2.8347637) q[0];
sx q[0];
rz(9.1260202) q[0];
rz(2.4755251) q[1];
sx q[1];
rz(-2.5112285) q[1];
sx q[1];
rz(-1.4730374) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10160343) q[0];
sx q[0];
rz(-2.371006) q[0];
sx q[0];
rz(-1.6165074) q[0];
rz(-pi) q[1];
rz(-0.017919964) q[2];
sx q[2];
rz(-1.8715053) q[2];
sx q[2];
rz(1.199388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.3949497) q[1];
sx q[1];
rz(-0.87374765) q[1];
sx q[1];
rz(-0.66956981) q[1];
x q[2];
rz(2.8657718) q[3];
sx q[3];
rz(-1.5032056) q[3];
sx q[3];
rz(1.3753613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95639688) q[2];
sx q[2];
rz(-2.7608725) q[2];
sx q[2];
rz(1.8308651) q[2];
rz(-3.0500566) q[3];
sx q[3];
rz(-2.5444578) q[3];
sx q[3];
rz(-2.6105647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700584) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(-0.48010215) q[0];
rz(-1.1473038) q[1];
sx q[1];
rz(-0.6310178) q[1];
sx q[1];
rz(0.70153418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7362721) q[0];
sx q[0];
rz(-0.87001409) q[0];
sx q[0];
rz(-2.2474656) q[0];
rz(0.417493) q[2];
sx q[2];
rz(-2.001363) q[2];
sx q[2];
rz(2.9232962) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50276668) q[1];
sx q[1];
rz(-2.354489) q[1];
sx q[1];
rz(-0.9197286) q[1];
x q[2];
rz(2.7357932) q[3];
sx q[3];
rz(-0.67849745) q[3];
sx q[3];
rz(-0.12756995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.39701617) q[2];
sx q[2];
rz(-1.1822367) q[2];
sx q[2];
rz(1.5783295) q[2];
rz(0.27966106) q[3];
sx q[3];
rz(-2.4249488) q[3];
sx q[3];
rz(0.40951148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79769832) q[0];
sx q[0];
rz(-0.41223031) q[0];
sx q[0];
rz(-1.718234) q[0];
rz(-0.1238981) q[1];
sx q[1];
rz(-0.84404498) q[1];
sx q[1];
rz(-1.557225) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6829872) q[0];
sx q[0];
rz(-1.2747972) q[0];
sx q[0];
rz(2.307555) q[0];
x q[1];
rz(-2.2940647) q[2];
sx q[2];
rz(-1.8516527) q[2];
sx q[2];
rz(2.6658863) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.023097087) q[1];
sx q[1];
rz(-1.6498936) q[1];
sx q[1];
rz(-1.3239856) q[1];
rz(-pi) q[2];
rz(-0.73488124) q[3];
sx q[3];
rz(-0.85414825) q[3];
sx q[3];
rz(-1.7295966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0763756) q[2];
sx q[2];
rz(-2.4880444) q[2];
sx q[2];
rz(0.93488133) q[2];
rz(-2.8377418) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(2.2936308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3728751) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(0.26707643) q[0];
rz(-3.0067387) q[1];
sx q[1];
rz(-0.57661533) q[1];
sx q[1];
rz(2.3042302) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42138153) q[0];
sx q[0];
rz(-2.1108339) q[0];
sx q[0];
rz(0.311894) q[0];
rz(-pi) q[1];
rz(2.1979273) q[2];
sx q[2];
rz(-0.84463464) q[2];
sx q[2];
rz(0.86792714) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3738651) q[1];
sx q[1];
rz(-0.30555913) q[1];
sx q[1];
rz(2.4453246) q[1];
rz(3.0947826) q[3];
sx q[3];
rz(-1.5093646) q[3];
sx q[3];
rz(0.84247473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59552938) q[2];
sx q[2];
rz(-2.9661621) q[2];
sx q[2];
rz(2.9626633) q[2];
rz(1.1692125) q[3];
sx q[3];
rz(-1.365265) q[3];
sx q[3];
rz(0.21807142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1201852) q[0];
sx q[0];
rz(-0.51161259) q[0];
sx q[0];
rz(0.87274337) q[0];
rz(-1.9841638) q[1];
sx q[1];
rz(-2.2667784) q[1];
sx q[1];
rz(0.016955888) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8411422) q[0];
sx q[0];
rz(-2.4068641) q[0];
sx q[0];
rz(-3.1233913) q[0];
rz(-pi) q[1];
rz(0.01466401) q[2];
sx q[2];
rz(-0.59978205) q[2];
sx q[2];
rz(1.5408915) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.03093623) q[1];
sx q[1];
rz(-0.82582322) q[1];
sx q[1];
rz(-0.40542545) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8722829) q[3];
sx q[3];
rz(-1.2896754) q[3];
sx q[3];
rz(-0.65015974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0852647) q[2];
sx q[2];
rz(-1.7076098) q[2];
sx q[2];
rz(-2.8223574) q[2];
rz(0.6790092) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(0.80176789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7769258) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(0.6426245) q[0];
rz(-1.369426) q[1];
sx q[1];
rz(-2.5182928) q[1];
sx q[1];
rz(-0.97698897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60212712) q[0];
sx q[0];
rz(-1.6940247) q[0];
sx q[0];
rz(1.6696403) q[0];
rz(-1.0426635) q[2];
sx q[2];
rz(-1.1234049) q[2];
sx q[2];
rz(1.3716979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.381666) q[1];
sx q[1];
rz(-1.5569512) q[1];
sx q[1];
rz(-2.2626876) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.053835458) q[3];
sx q[3];
rz(-2.1559058) q[3];
sx q[3];
rz(-1.9991635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6151108) q[2];
sx q[2];
rz(-1.7340163) q[2];
sx q[2];
rz(-1.2283481) q[2];
rz(2.2743716) q[3];
sx q[3];
rz(-2.7286178) q[3];
sx q[3];
rz(-0.76798463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41384554) q[0];
sx q[0];
rz(-2.9242046) q[0];
sx q[0];
rz(-2.9687498) q[0];
rz(-0.8051644) q[1];
sx q[1];
rz(-2.5925437) q[1];
sx q[1];
rz(-0.13596143) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8589218) q[0];
sx q[0];
rz(-1.4513119) q[0];
sx q[0];
rz(-0.8651328) q[0];
rz(-pi) q[1];
rz(0.82201652) q[2];
sx q[2];
rz(-1.4944585) q[2];
sx q[2];
rz(-0.99886307) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0522642) q[1];
sx q[1];
rz(-1.58984) q[1];
sx q[1];
rz(-0.11472265) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37660677) q[3];
sx q[3];
rz(-2.7889502) q[3];
sx q[3];
rz(2.6654748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92676306) q[2];
sx q[2];
rz(-1.5952933) q[2];
sx q[2];
rz(0.1845486) q[2];
rz(-0.18937011) q[3];
sx q[3];
rz(-2.6034077) q[3];
sx q[3];
rz(-2.2034933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1592584) q[0];
sx q[0];
rz(-0.33894798) q[0];
sx q[0];
rz(-2.2151997) q[0];
rz(3.1023846) q[1];
sx q[1];
rz(-2.4640633) q[1];
sx q[1];
rz(2.1047986) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2873357) q[0];
sx q[0];
rz(-1.8651891) q[0];
sx q[0];
rz(0.18369412) q[0];
rz(-2.762297) q[2];
sx q[2];
rz(-2.6808028) q[2];
sx q[2];
rz(-0.57578218) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29940393) q[1];
sx q[1];
rz(-1.328598) q[1];
sx q[1];
rz(-2.7224225) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7091712) q[3];
sx q[3];
rz(-2.7161971) q[3];
sx q[3];
rz(-2.2063125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2490273) q[2];
sx q[2];
rz(-1.1250863) q[2];
sx q[2];
rz(0.79891515) q[2];
rz(-0.51698452) q[3];
sx q[3];
rz(-2.7345782) q[3];
sx q[3];
rz(-0.5504722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.35140458) q[0];
sx q[0];
rz(-3.1365972) q[0];
sx q[0];
rz(-1.6222401) q[0];
rz(1.5478569) q[1];
sx q[1];
rz(-2.3212815) q[1];
sx q[1];
rz(-0.69040745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0102745) q[0];
sx q[0];
rz(-1.2206435) q[0];
sx q[0];
rz(1.8084333) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0418441) q[2];
sx q[2];
rz(-2.3398551) q[2];
sx q[2];
rz(-0.24101098) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71296299) q[1];
sx q[1];
rz(-2.3979514) q[1];
sx q[1];
rz(2.3548467) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9966808) q[3];
sx q[3];
rz(-2.1050801) q[3];
sx q[3];
rz(-0.89569672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7030299) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(-2.1717066) q[2];
rz(0.25740933) q[3];
sx q[3];
rz(-2.9195547) q[3];
sx q[3];
rz(-1.7820057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6064706) q[0];
sx q[0];
rz(-0.73000014) q[0];
sx q[0];
rz(-0.085513376) q[0];
rz(2.8657148) q[1];
sx q[1];
rz(-2.52067) q[1];
sx q[1];
rz(1.9524908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71578854) q[0];
sx q[0];
rz(-0.5606519) q[0];
sx q[0];
rz(0.7348357) q[0];
x q[1];
rz(0.3100119) q[2];
sx q[2];
rz(-1.5622721) q[2];
sx q[2];
rz(-0.046176813) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0240636) q[1];
sx q[1];
rz(-1.4520565) q[1];
sx q[1];
rz(-1.9734782) q[1];
x q[2];
rz(-2.7932634) q[3];
sx q[3];
rz(-1.4831721) q[3];
sx q[3];
rz(-1.5867811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6659866) q[2];
sx q[2];
rz(-0.87916547) q[2];
sx q[2];
rz(0.48412588) q[2];
rz(-0.13949805) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(-2.9805396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465268) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(1.4576661) q[1];
sx q[1];
rz(-1.6438345) q[1];
sx q[1];
rz(1.8297292) q[1];
rz(-0.64907907) q[2];
sx q[2];
rz(-1.4076283) q[2];
sx q[2];
rz(0.84244737) q[2];
rz(-3.1112989) q[3];
sx q[3];
rz(-0.75728635) q[3];
sx q[3];
rz(-1.5273619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
