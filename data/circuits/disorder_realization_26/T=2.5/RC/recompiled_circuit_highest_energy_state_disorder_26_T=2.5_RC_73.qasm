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
rz(0.8910203) q[0];
sx q[0];
rz(-1.2863337) q[0];
sx q[0];
rz(0.88598716) q[0];
rz(-2.6180144) q[1];
sx q[1];
rz(-0.55966592) q[1];
sx q[1];
rz(-1.8212512) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26368982) q[0];
sx q[0];
rz(-0.56803507) q[0];
sx q[0];
rz(-1.2391511) q[0];
x q[1];
rz(2.7971326) q[2];
sx q[2];
rz(-0.80756029) q[2];
sx q[2];
rz(-1.5429614) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3229827) q[1];
sx q[1];
rz(-0.56639379) q[1];
sx q[1];
rz(1.4569805) q[1];
rz(2.376972) q[3];
sx q[3];
rz(-0.9249827) q[3];
sx q[3];
rz(-0.82987204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77871275) q[2];
sx q[2];
rz(-1.820463) q[2];
sx q[2];
rz(0.92817456) q[2];
rz(-1.5859531) q[3];
sx q[3];
rz(-1.0471683) q[3];
sx q[3];
rz(1.1699404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8856119) q[0];
sx q[0];
rz(-2.8035127) q[0];
sx q[0];
rz(2.0804491) q[0];
rz(-1.9288918) q[1];
sx q[1];
rz(-2.3586528) q[1];
sx q[1];
rz(-1.8276259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3201225) q[0];
sx q[0];
rz(-2.263453) q[0];
sx q[0];
rz(-2.0351324) q[0];
rz(-pi) q[1];
rz(2.7065837) q[2];
sx q[2];
rz(-0.97112964) q[2];
sx q[2];
rz(0.51520447) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.48294327) q[1];
sx q[1];
rz(-1.4023781) q[1];
sx q[1];
rz(-1.0337109) q[1];
x q[2];
rz(-2.1589018) q[3];
sx q[3];
rz(-0.81702166) q[3];
sx q[3];
rz(-0.52023653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5709915) q[2];
sx q[2];
rz(-2.1676895) q[2];
sx q[2];
rz(2.2368597) q[2];
rz(1.5500801) q[3];
sx q[3];
rz(-1.855987) q[3];
sx q[3];
rz(-1.8486283) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631113) q[0];
sx q[0];
rz(-3.0516629) q[0];
sx q[0];
rz(2.2233546) q[0];
rz(-2.4507484) q[1];
sx q[1];
rz(-1.7450688) q[1];
sx q[1];
rz(-2.3450559) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18609234) q[0];
sx q[0];
rz(-1.4791795) q[0];
sx q[0];
rz(1.7797526) q[0];
rz(2.1941625) q[2];
sx q[2];
rz(-0.073270144) q[2];
sx q[2];
rz(0.87033349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0738363) q[1];
sx q[1];
rz(-2.2646565) q[1];
sx q[1];
rz(0.14090726) q[1];
rz(1.1311319) q[3];
sx q[3];
rz(-0.70207233) q[3];
sx q[3];
rz(2.0560101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1306661) q[2];
sx q[2];
rz(-2.1911759) q[2];
sx q[2];
rz(1.8939023) q[2];
rz(0.55825663) q[3];
sx q[3];
rz(-1.4562675) q[3];
sx q[3];
rz(0.021350967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.438544) q[0];
sx q[0];
rz(-0.049228638) q[0];
sx q[0];
rz(-2.4141648) q[0];
rz(2.9797331) q[1];
sx q[1];
rz(-1.1715803) q[1];
sx q[1];
rz(-1.9836099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25140554) q[0];
sx q[0];
rz(-1.5047538) q[0];
sx q[0];
rz(1.6679881) q[0];
rz(-pi) q[1];
rz(1.7100934) q[2];
sx q[2];
rz(-1.5717531) q[2];
sx q[2];
rz(0.19569163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47745142) q[1];
sx q[1];
rz(-0.22813561) q[1];
sx q[1];
rz(2.1456108) q[1];
x q[2];
rz(-1.653966) q[3];
sx q[3];
rz(-1.0462049) q[3];
sx q[3];
rz(0.34357854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27266476) q[2];
sx q[2];
rz(-1.9395892) q[2];
sx q[2];
rz(-1.7499917) q[2];
rz(2.9113655) q[3];
sx q[3];
rz(-1.1743436) q[3];
sx q[3];
rz(0.19190425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4153083) q[0];
sx q[0];
rz(-3.0400161) q[0];
sx q[0];
rz(-2.0368077) q[0];
rz(3.0305908) q[1];
sx q[1];
rz(-1.1155201) q[1];
sx q[1];
rz(1.312779) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5376624) q[0];
sx q[0];
rz(-0.27034187) q[0];
sx q[0];
rz(2.2946847) q[0];
rz(1.6592624) q[2];
sx q[2];
rz(-1.0087722) q[2];
sx q[2];
rz(-2.1254181) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.032636383) q[1];
sx q[1];
rz(-1.3707038) q[1];
sx q[1];
rz(1.97831) q[1];
rz(-pi) q[2];
rz(2.2897374) q[3];
sx q[3];
rz(-1.2491711) q[3];
sx q[3];
rz(0.77675288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.50292) q[2];
sx q[2];
rz(-1.0588249) q[2];
sx q[2];
rz(1.0271094) q[2];
rz(2.1549759) q[3];
sx q[3];
rz(-2.8592181) q[3];
sx q[3];
rz(-1.4237039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33243772) q[0];
sx q[0];
rz(-1.2323392) q[0];
sx q[0];
rz(-0.80068457) q[0];
rz(-2.370131) q[1];
sx q[1];
rz(-1.0779251) q[1];
sx q[1];
rz(1.3763743) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0906935) q[0];
sx q[0];
rz(-0.92779033) q[0];
sx q[0];
rz(1.0674632) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.034997392) q[2];
sx q[2];
rz(-2.0366324) q[2];
sx q[2];
rz(-2.310487) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.47471646) q[1];
sx q[1];
rz(-2.2762269) q[1];
sx q[1];
rz(1.553346) q[1];
rz(-2.3153051) q[3];
sx q[3];
rz(-1.0184231) q[3];
sx q[3];
rz(2.1635522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8419522) q[2];
sx q[2];
rz(-1.9349293) q[2];
sx q[2];
rz(0.017814962) q[2];
rz(-0.31339112) q[3];
sx q[3];
rz(-1.0217977) q[3];
sx q[3];
rz(-0.71948403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.3571091) q[0];
sx q[0];
rz(-1.5743558) q[0];
sx q[0];
rz(-0.16726476) q[0];
rz(0.74370614) q[1];
sx q[1];
rz(-2.2552762) q[1];
sx q[1];
rz(3.0053265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1377099) q[0];
sx q[0];
rz(-1.2943177) q[0];
sx q[0];
rz(2.2599758) q[0];
rz(-pi) q[1];
rz(1.0616779) q[2];
sx q[2];
rz(-2.7973242) q[2];
sx q[2];
rz(1.6580251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2124575) q[1];
sx q[1];
rz(-1.2127185) q[1];
sx q[1];
rz(1.6422988) q[1];
rz(-pi) q[2];
rz(-1.6266009) q[3];
sx q[3];
rz(-1.2254834) q[3];
sx q[3];
rz(-1.8453516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3027489) q[2];
sx q[2];
rz(-2.9911797) q[2];
sx q[2];
rz(-2.0602843) q[2];
rz(2.9980764) q[3];
sx q[3];
rz(-1.2986526) q[3];
sx q[3];
rz(-1.6941841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3926587) q[0];
sx q[0];
rz(-1.8086139) q[0];
sx q[0];
rz(-0.64312154) q[0];
rz(2.2195623) q[1];
sx q[1];
rz(-2.8755201) q[1];
sx q[1];
rz(1.6304852) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6623343) q[0];
sx q[0];
rz(-1.3915724) q[0];
sx q[0];
rz(-1.2249806) q[0];
x q[1];
rz(3.0273075) q[2];
sx q[2];
rz(-1.8488036) q[2];
sx q[2];
rz(2.4604527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1280316) q[1];
sx q[1];
rz(-1.7674539) q[1];
sx q[1];
rz(0.0081045919) q[1];
rz(-0.51120843) q[3];
sx q[3];
rz(-1.852559) q[3];
sx q[3];
rz(0.35387999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1954605) q[2];
sx q[2];
rz(-1.5798502) q[2];
sx q[2];
rz(1.3004318) q[2];
rz(2.2293633) q[3];
sx q[3];
rz(-1.8650841) q[3];
sx q[3];
rz(0.31360489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6533971) q[0];
sx q[0];
rz(-0.47976872) q[0];
sx q[0];
rz(2.524014) q[0];
rz(-0.081534475) q[1];
sx q[1];
rz(-2.640994) q[1];
sx q[1];
rz(0.68731442) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36022913) q[0];
sx q[0];
rz(-1.5232674) q[0];
sx q[0];
rz(1.8076872) q[0];
x q[1];
rz(0.26726802) q[2];
sx q[2];
rz(-1.6785926) q[2];
sx q[2];
rz(1.9330618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1256492) q[1];
sx q[1];
rz(-1.9852891) q[1];
sx q[1];
rz(-1.260956) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1140442) q[3];
sx q[3];
rz(-1.5294269) q[3];
sx q[3];
rz(-0.70391212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0240137) q[2];
sx q[2];
rz(-1.1774096) q[2];
sx q[2];
rz(-0.071361072) q[2];
rz(1.2355545) q[3];
sx q[3];
rz(-0.33696285) q[3];
sx q[3];
rz(0.56628847) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5736893) q[0];
sx q[0];
rz(-2.8901143) q[0];
sx q[0];
rz(0.13791826) q[0];
rz(2.5714286) q[1];
sx q[1];
rz(-1.0823931) q[1];
sx q[1];
rz(-0.015448419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.954166) q[0];
sx q[0];
rz(-1.2492234) q[0];
sx q[0];
rz(-2.2963752) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49039109) q[2];
sx q[2];
rz(-0.82329673) q[2];
sx q[2];
rz(3.1237912) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9271157) q[1];
sx q[1];
rz(-2.0600015) q[1];
sx q[1];
rz(1.3838883) q[1];
rz(-pi) q[2];
rz(1.3872765) q[3];
sx q[3];
rz(-1.6257111) q[3];
sx q[3];
rz(-0.6764937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7323759) q[2];
sx q[2];
rz(-0.78745431) q[2];
sx q[2];
rz(2.709205) q[2];
rz(2.2435097) q[3];
sx q[3];
rz(-2.2663074) q[3];
sx q[3];
rz(-1.5842277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4846004) q[0];
sx q[0];
rz(-2.0130172) q[0];
sx q[0];
rz(0.76242557) q[0];
rz(-2.5905329) q[1];
sx q[1];
rz(-1.8012128) q[1];
sx q[1];
rz(2.6691379) q[1];
rz(-0.1189143) q[2];
sx q[2];
rz(-2.7758895) q[2];
sx q[2];
rz(2.9056673) q[2];
rz(-2.1848706) q[3];
sx q[3];
rz(-1.9695558) q[3];
sx q[3];
rz(1.1576681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
