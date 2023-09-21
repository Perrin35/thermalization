OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(-2.9876246) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(-2.8032803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96903893) q[0];
sx q[0];
rz(-1.8882897) q[0];
sx q[0];
rz(-2.88455) q[0];
rz(-1.1233653) q[2];
sx q[2];
rz(-2.4123203) q[2];
sx q[2];
rz(-2.4016618) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0677883) q[1];
sx q[1];
rz(-1.6238302) q[1];
sx q[1];
rz(-2.8552613) q[1];
rz(3.1258718) q[3];
sx q[3];
rz(-2.0831046) q[3];
sx q[3];
rz(-2.6920464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-0.34037408) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(-3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3409815) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(-0.064963438) q[0];
rz(-0.57463542) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(1.8992791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9984098) q[0];
sx q[0];
rz(-1.3230723) q[0];
sx q[0];
rz(-1.7032911) q[0];
x q[1];
rz(-1.3556446) q[2];
sx q[2];
rz(-0.64385507) q[2];
sx q[2];
rz(-1.3607963) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4878792) q[1];
sx q[1];
rz(-2.5839845) q[1];
sx q[1];
rz(-2.8370268) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70988016) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(-0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3339281) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(2.5578965) q[2];
rz(-2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7211001) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(-0.72845355) q[0];
rz(1.6473673) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(1.0167936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66660488) q[0];
sx q[0];
rz(-3.0525065) q[0];
sx q[0];
rz(2.7276917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.617308) q[2];
sx q[2];
rz(-1.1384083) q[2];
sx q[2];
rz(1.0334894) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96117561) q[1];
sx q[1];
rz(-2.3824586) q[1];
sx q[1];
rz(-1.3707861) q[1];
x q[2];
rz(-1.4025027) q[3];
sx q[3];
rz(-2.8561391) q[3];
sx q[3];
rz(1.9608378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8016522) q[2];
sx q[2];
rz(-1.5321956) q[2];
sx q[2];
rz(-2.0920848) q[2];
rz(2.5028051) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(-1.9558186) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8124354) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.4720434) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(0.24681117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1082552) q[0];
sx q[0];
rz(-0.70972432) q[0];
sx q[0];
rz(1.8002585) q[0];
x q[1];
rz(-0.47841448) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(1.071655) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3226763) q[1];
sx q[1];
rz(-2.7671742) q[1];
sx q[1];
rz(-0.14426343) q[1];
rz(0.95440063) q[3];
sx q[3];
rz(-1.8064926) q[3];
sx q[3];
rz(0.97782048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4776769) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33070579) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(-1.8141618) q[0];
rz(-1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(2.8932103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029593) q[0];
sx q[0];
rz(-2.6002433) q[0];
sx q[0];
rz(0.066141733) q[0];
x q[1];
rz(0.70154538) q[2];
sx q[2];
rz(-2.8994459) q[2];
sx q[2];
rz(2.0602351) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0866962) q[1];
sx q[1];
rz(-2.5270562) q[1];
sx q[1];
rz(2.8413248) q[1];
x q[2];
rz(-1.3316657) q[3];
sx q[3];
rz(-0.93591792) q[3];
sx q[3];
rz(-2.2374416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7053232) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(0.38875368) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-2.6388772) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0080863) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.3457993) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(0.1246917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69783083) q[0];
sx q[0];
rz(-2.9450581) q[0];
sx q[0];
rz(0.4561119) q[0];
rz(-2.5325534) q[2];
sx q[2];
rz(-1.3281203) q[2];
sx q[2];
rz(2.2720624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5387419) q[1];
sx q[1];
rz(-2.3067143) q[1];
sx q[1];
rz(-1.49453) q[1];
rz(-pi) q[2];
rz(2.9051404) q[3];
sx q[3];
rz(-1.3580305) q[3];
sx q[3];
rz(-2.3355683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51320118) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(-0.61895269) q[2];
rz(-1.0533054) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(0.9296023) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182794) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(2.0986309) q[0];
rz(-2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(-2.0708864) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8541504) q[0];
sx q[0];
rz(-1.5713912) q[0];
sx q[0];
rz(2.6575412) q[0];
rz(2.7878739) q[2];
sx q[2];
rz(-0.36834799) q[2];
sx q[2];
rz(3.0596717) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5412022) q[1];
sx q[1];
rz(-1.6704847) q[1];
sx q[1];
rz(0.77460918) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.058810874) q[3];
sx q[3];
rz(-2.8088514) q[3];
sx q[3];
rz(-2.9174093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-0.32361844) q[2];
rz(-0.98179022) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535646) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(-1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(0.94747296) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5763801) q[0];
sx q[0];
rz(-2.1064261) q[0];
sx q[0];
rz(-0.11760786) q[0];
rz(1.7296373) q[2];
sx q[2];
rz(-1.7454595) q[2];
sx q[2];
rz(0.58089248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35093388) q[1];
sx q[1];
rz(-0.62218636) q[1];
sx q[1];
rz(-1.15637) q[1];
rz(2.794907) q[3];
sx q[3];
rz(-0.83804916) q[3];
sx q[3];
rz(0.29688641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(2.6718111) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(-0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(-0.20283094) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79265362) q[0];
sx q[0];
rz(-1.5447504) q[0];
sx q[0];
rz(-3.1025725) q[0];
rz(2.3896396) q[2];
sx q[2];
rz(-2.8402036) q[2];
sx q[2];
rz(1.6981268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54535941) q[1];
sx q[1];
rz(-1.4701478) q[1];
sx q[1];
rz(-0.049493162) q[1];
rz(-pi) q[2];
rz(-2.6796954) q[3];
sx q[3];
rz(-1.4326722) q[3];
sx q[3];
rz(-1.671333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7245076) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-2.1450796) q[2];
rz(0.35342446) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(-0.18173519) q[0];
rz(-3.0985447) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(-0.28082401) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1998394) q[0];
sx q[0];
rz(-1.117525) q[0];
sx q[0];
rz(2.0487294) q[0];
rz(-pi) q[1];
rz(-1.3781204) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(-0.10373058) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2320497) q[1];
sx q[1];
rz(-0.068040158) q[1];
sx q[1];
rz(1.0104936) q[1];
rz(0.19631581) q[3];
sx q[3];
rz(-1.2442949) q[3];
sx q[3];
rz(1.452009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3165555) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(0.62310702) q[2];
rz(2.1394219) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(-0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(2.2676246) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-2.5429824) q[2];
sx q[2];
rz(-2.241588) q[2];
sx q[2];
rz(-0.66551756) q[2];
rz(-1.9813886) q[3];
sx q[3];
rz(-1.900233) q[3];
sx q[3];
rz(-1.486447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];