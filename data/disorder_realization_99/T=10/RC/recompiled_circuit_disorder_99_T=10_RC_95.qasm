OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0302304) q[0];
sx q[0];
rz(5.6279866) q[0];
sx q[0];
rz(10.401166) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3800669) q[0];
sx q[0];
rz(-0.67136429) q[0];
sx q[0];
rz(-1.5321561) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36926271) q[2];
sx q[2];
rz(-2.5662176) q[2];
sx q[2];
rz(2.0865666) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9470889) q[1];
sx q[1];
rz(-2.6852486) q[1];
sx q[1];
rz(1.7366921) q[1];
rz(-2.3402432) q[3];
sx q[3];
rz(-2.1270463) q[3];
sx q[3];
rz(0.43503161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4684326) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(-0.012714816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760261) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-0.20543988) q[1];
sx q[1];
rz(-0.36546779) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.394709) q[0];
sx q[0];
rz(-0.53046662) q[0];
sx q[0];
rz(2.6283426) q[0];
x q[1];
rz(-2.2720049) q[2];
sx q[2];
rz(-0.80539942) q[2];
sx q[2];
rz(2.373901) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47479113) q[1];
sx q[1];
rz(-2.431369) q[1];
sx q[1];
rz(2.6067961) q[1];
rz(2.6037381) q[3];
sx q[3];
rz(-2.9085277) q[3];
sx q[3];
rz(-0.92249289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83313292) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66115528) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(0.61082947) q[0];
rz(1.4340596) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(-2.6909713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47942802) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(1.9730294) q[0];
x q[1];
rz(0.30655105) q[2];
sx q[2];
rz(-1.2129285) q[2];
sx q[2];
rz(-1.7265665) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9420942) q[1];
sx q[1];
rz(-2.4220785) q[1];
sx q[1];
rz(0.47902963) q[1];
rz(-pi) q[2];
rz(1.9335453) q[3];
sx q[3];
rz(-0.72788531) q[3];
sx q[3];
rz(-0.7286275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11015636) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(2.2971161) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-2.0003831) q[0];
rz(-1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(-0.30532125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13320623) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(-0.58831373) q[0];
rz(1.1496454) q[2];
sx q[2];
rz(-1.1425758) q[2];
sx q[2];
rz(2.2575833) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1034643) q[1];
sx q[1];
rz(-2.3410019) q[1];
sx q[1];
rz(2.994328) q[1];
rz(2.6257315) q[3];
sx q[3];
rz(-2.0277884) q[3];
sx q[3];
rz(-2.7716694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9332463) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(2.7944881) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(-0.53891671) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(2.7358823) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1203128) q[0];
sx q[0];
rz(-2.4097674) q[0];
sx q[0];
rz(0.41723199) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0056555) q[2];
sx q[2];
rz(-1.3456151) q[2];
sx q[2];
rz(-1.2852247) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5647445) q[1];
sx q[1];
rz(-1.2291359) q[1];
sx q[1];
rz(0.26775743) q[1];
x q[2];
rz(-1.8040854) q[3];
sx q[3];
rz(-0.81987112) q[3];
sx q[3];
rz(-1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47200176) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(-1.4542788) q[2];
rz(-0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-2.7344761) q[0];
rz(0.85917568) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28535298) q[0];
sx q[0];
rz(-1.0278773) q[0];
sx q[0];
rz(3.1154484) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98951927) q[2];
sx q[2];
rz(-2.1534854) q[2];
sx q[2];
rz(-1.4942102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0504005) q[1];
sx q[1];
rz(-0.7351774) q[1];
sx q[1];
rz(-1.5749732) q[1];
x q[2];
rz(-2.7139211) q[3];
sx q[3];
rz(-1.0276405) q[3];
sx q[3];
rz(0.021051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(1.9334531) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5642501) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(-1.1279001) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.8022585) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3198277) q[0];
sx q[0];
rz(-1.9263096) q[0];
sx q[0];
rz(-2.6000644) q[0];
rz(-2.5008051) q[2];
sx q[2];
rz(-1.0377585) q[2];
sx q[2];
rz(1.2696881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5994508) q[1];
sx q[1];
rz(-2.8891924) q[1];
sx q[1];
rz(-0.42599328) q[1];
x q[2];
rz(-3.1197458) q[3];
sx q[3];
rz(-2.1751746) q[3];
sx q[3];
rz(1.5165839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82905519) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(-2.6230295) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(-1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.1475875) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(2.5491319) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5151278) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(-1.5149087) q[0];
rz(-pi) q[1];
rz(-0.7166491) q[2];
sx q[2];
rz(-0.71084329) q[2];
sx q[2];
rz(-0.58472842) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8471624) q[1];
sx q[1];
rz(-1.2420328) q[1];
sx q[1];
rz(-0.86818236) q[1];
rz(1.1660277) q[3];
sx q[3];
rz(-1.3169857) q[3];
sx q[3];
rz(0.39241957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56117326) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(-2.6994761) q[2];
rz(2.2413065) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(-3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(1.3046718) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-2.5794199) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32414831) q[0];
sx q[0];
rz(-1.9262909) q[0];
sx q[0];
rz(-0.72974156) q[0];
rz(2.1015342) q[2];
sx q[2];
rz(-1.4897523) q[2];
sx q[2];
rz(-2.9329252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1162029) q[1];
sx q[1];
rz(-0.89239489) q[1];
sx q[1];
rz(-2.5118026) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1702873) q[3];
sx q[3];
rz(-1.1937965) q[3];
sx q[3];
rz(-1.4124944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(0.98158681) q[2];
rz(3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-0.04709588) q[0];
rz(-1.3866407) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(-0.047853619) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82057792) q[0];
sx q[0];
rz(-0.65549675) q[0];
sx q[0];
rz(-1.5232248) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1282975) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(-2.0890582) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25803265) q[1];
sx q[1];
rz(-1.1415392) q[1];
sx q[1];
rz(-0.92991035) q[1];
rz(2.738897) q[3];
sx q[3];
rz(-1.5201357) q[3];
sx q[3];
rz(2.5718873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(0.93635526) q[2];
rz(0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8235648) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(0.035632523) q[2];
sx q[2];
rz(-2.590066) q[2];
sx q[2];
rz(-1.8077015) q[2];
rz(1.1669284) q[3];
sx q[3];
rz(-1.231791) q[3];
sx q[3];
rz(0.024699208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
