OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(0.18145951) q[0];
rz(2.060086) q[1];
sx q[1];
rz(5.6097538) q[1];
sx q[1];
rz(14.654832) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5232789) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(2.0413415) q[0];
rz(-pi) q[1];
rz(-2.7396169) q[2];
sx q[2];
rz(-0.66265124) q[2];
sx q[2];
rz(2.0768009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3095113) q[1];
sx q[1];
rz(-0.81334269) q[1];
sx q[1];
rz(1.3002212) q[1];
rz(-pi) q[2];
rz(1.1374723) q[3];
sx q[3];
rz(-2.7702799) q[3];
sx q[3];
rz(0.565688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(-1.0129499) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-0.071454123) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.6275303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1925416) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(0.87829725) q[0];
rz(-1.6018277) q[2];
sx q[2];
rz(-0.73756644) q[2];
sx q[2];
rz(2.8013128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82799998) q[1];
sx q[1];
rz(-0.26121059) q[1];
sx q[1];
rz(2.2952495) q[1];
x q[2];
rz(1.1176923) q[3];
sx q[3];
rz(-3.0278904) q[3];
sx q[3];
rz(0.5773471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4864768) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(-0.97186175) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(-2.8308716) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(-1.4583189) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70143914) q[0];
sx q[0];
rz(-2.9868786) q[0];
sx q[0];
rz(-1.8811149) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4086191) q[2];
sx q[2];
rz(-0.79967116) q[2];
sx q[2];
rz(3.1250931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2787784) q[1];
sx q[1];
rz(-1.9796951) q[1];
sx q[1];
rz(-1.6294953) q[1];
rz(-pi) q[2];
rz(-0.31483105) q[3];
sx q[3];
rz(-1.587084) q[3];
sx q[3];
rz(0.44486886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(0.67908755) q[2];
rz(1.7689765) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2963592) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.4845928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35614466) q[0];
sx q[0];
rz(-1.5404285) q[0];
sx q[0];
rz(-1.5201475) q[0];
rz(-pi) q[1];
rz(1.053327) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(0.53211624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1160108) q[1];
sx q[1];
rz(-0.75041795) q[1];
sx q[1];
rz(0.90390272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4184166) q[3];
sx q[3];
rz(-1.8225065) q[3];
sx q[3];
rz(0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(1.0162639) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6338585) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(-0.39988363) q[0];
rz(1.9790861) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(-2.9679325) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6337316) q[0];
sx q[0];
rz(-1.5426239) q[0];
sx q[0];
rz(-0.072576056) q[0];
rz(-pi) q[1];
rz(-2.747614) q[2];
sx q[2];
rz(-1.3756961) q[2];
sx q[2];
rz(-0.56001284) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34827161) q[1];
sx q[1];
rz(-1.4655359) q[1];
sx q[1];
rz(3.1131016) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3624304) q[3];
sx q[3];
rz(-1.2456129) q[3];
sx q[3];
rz(1.9735826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48012039) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(2.373467) q[2];
rz(0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(-0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0356692) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(1.6712028) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(1.1434198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238732) q[0];
sx q[0];
rz(-1.4871162) q[0];
sx q[0];
rz(0.15079817) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46164718) q[2];
sx q[2];
rz(-2.4857773) q[2];
sx q[2];
rz(2.0222424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.45149849) q[1];
sx q[1];
rz(-1.203981) q[1];
sx q[1];
rz(-2.9551798) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8361736) q[3];
sx q[3];
rz(-1.050204) q[3];
sx q[3];
rz(0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(-1.6879843) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.24916515) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221508) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(2.4123736) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-1.1332606) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6183375) q[0];
sx q[0];
rz(-1.6132857) q[0];
sx q[0];
rz(2.8190814) q[0];
rz(-pi) q[1];
rz(0.015958162) q[2];
sx q[2];
rz(-0.71231132) q[2];
sx q[2];
rz(2.6238837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3164191) q[1];
sx q[1];
rz(-0.9269886) q[1];
sx q[1];
rz(-0.34367798) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4134737) q[3];
sx q[3];
rz(-1.8295349) q[3];
sx q[3];
rz(1.7595921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7363654) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(-0.48842946) q[2];
rz(1.8296261) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(-0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768123) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(-0.07117614) q[0];
rz(-0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.9326928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29584822) q[0];
sx q[0];
rz(-0.83202067) q[0];
sx q[0];
rz(2.8315298) q[0];
rz(2.2279943) q[2];
sx q[2];
rz(-2.845394) q[2];
sx q[2];
rz(-1.9343455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88747957) q[1];
sx q[1];
rz(-1.5548318) q[1];
sx q[1];
rz(1.70114) q[1];
rz(-1.9367996) q[3];
sx q[3];
rz(-0.98099785) q[3];
sx q[3];
rz(-2.3596711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9341087) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(-0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97380012) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(0.31162509) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.5100381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2841543) q[0];
sx q[0];
rz(-2.7630002) q[0];
sx q[0];
rz(2.5749102) q[0];
rz(-pi) q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(2.0575112) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9540625) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(0.33893434) q[1];
rz(2.4771677) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2733549) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(-2.3256425) q[2];
rz(-2.6319035) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508535) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(-2.5157805) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(-0.65840107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9068245) q[0];
sx q[0];
rz(-1.0848197) q[0];
sx q[0];
rz(-2.4776239) q[0];
rz(-pi) q[1];
rz(-0.41215956) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(1.5199496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8734332) q[1];
sx q[1];
rz(-2.6323316) q[1];
sx q[1];
rz(1.9164273) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9202616) q[3];
sx q[3];
rz(-1.06171) q[3];
sx q[3];
rz(1.6125319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4663503) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(-2.8038483) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.52453775) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(0.94454371) q[2];
sx q[2];
rz(-1.4251475) q[2];
sx q[2];
rz(0.057694358) q[2];
rz(-3.0425439) q[3];
sx q[3];
rz(-0.5746114) q[3];
sx q[3];
rz(-0.050669908) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];