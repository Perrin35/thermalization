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
rz(1.5647178) q[0];
sx q[0];
rz(-1.325542) q[0];
sx q[0];
rz(2.5817885) q[0];
rz(-1.7436104) q[1];
sx q[1];
rz(1.8140732) q[1];
sx q[1];
rz(10.811364) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35980269) q[0];
sx q[0];
rz(-1.3179147) q[0];
sx q[0];
rz(2.1232312) q[0];
x q[1];
rz(-0.29560372) q[2];
sx q[2];
rz(-1.6787375) q[2];
sx q[2];
rz(-2.2578277) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0530449) q[1];
sx q[1];
rz(-1.7783718) q[1];
sx q[1];
rz(-1.6563841) q[1];
rz(-pi) q[2];
rz(1.54744) q[3];
sx q[3];
rz(-0.29072116) q[3];
sx q[3];
rz(-0.61222044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5379415) q[2];
sx q[2];
rz(-0.62749046) q[2];
sx q[2];
rz(0.44169912) q[2];
rz(-0.80080992) q[3];
sx q[3];
rz(-1.2468612) q[3];
sx q[3];
rz(-2.7504564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(1.2342728) q[0];
sx q[0];
rz(-1.3751625) q[0];
sx q[0];
rz(2.8148839) q[0];
rz(1.3455343) q[1];
sx q[1];
rz(-1.5252472) q[1];
sx q[1];
rz(-2.2950642) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.617793) q[0];
sx q[0];
rz(-1.9547858) q[0];
sx q[0];
rz(0.11597534) q[0];
rz(-pi) q[1];
rz(-0.2791762) q[2];
sx q[2];
rz(-0.90718944) q[2];
sx q[2];
rz(0.070726591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1364899) q[1];
sx q[1];
rz(-2.0999036) q[1];
sx q[1];
rz(0.012261475) q[1];
rz(-pi) q[2];
rz(-1.2760824) q[3];
sx q[3];
rz(-0.8407774) q[3];
sx q[3];
rz(1.7129912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.232406) q[2];
sx q[2];
rz(-1.0973944) q[2];
sx q[2];
rz(0.48420134) q[2];
rz(1.3062612) q[3];
sx q[3];
rz(-2.5808915) q[3];
sx q[3];
rz(-1.1854712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4803798) q[0];
sx q[0];
rz(-0.73037195) q[0];
sx q[0];
rz(0.54006201) q[0];
rz(1.9909319) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(0.59404341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4503419) q[0];
sx q[0];
rz(-1.6800386) q[0];
sx q[0];
rz(1.3997549) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8833022) q[2];
sx q[2];
rz(-2.2021026) q[2];
sx q[2];
rz(-2.3986651) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3004527) q[1];
sx q[1];
rz(-2.3256635) q[1];
sx q[1];
rz(-3.0197179) q[1];
rz(1.0888799) q[3];
sx q[3];
rz(-1.6301117) q[3];
sx q[3];
rz(2.7869567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8192886) q[2];
sx q[2];
rz(-2.1257336) q[2];
sx q[2];
rz(-2.4508396) q[2];
rz(-0.47762075) q[3];
sx q[3];
rz(-2.728929) q[3];
sx q[3];
rz(-0.42948183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531042) q[0];
sx q[0];
rz(-0.34585837) q[0];
sx q[0];
rz(-2.389287) q[0];
rz(0.90657702) q[1];
sx q[1];
rz(-2.7598359) q[1];
sx q[1];
rz(0.22901542) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54628263) q[0];
sx q[0];
rz(-2.7685099) q[0];
sx q[0];
rz(0.65600106) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7976409) q[2];
sx q[2];
rz(-2.0984435) q[2];
sx q[2];
rz(-0.91560748) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6743857) q[1];
sx q[1];
rz(-1.8388646) q[1];
sx q[1];
rz(1.3246791) q[1];
rz(-pi) q[2];
rz(-2.6182272) q[3];
sx q[3];
rz(-0.65306758) q[3];
sx q[3];
rz(-3.0619597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0166246) q[2];
sx q[2];
rz(-0.11398537) q[2];
sx q[2];
rz(-1.0329186) q[2];
rz(-1.0659263) q[3];
sx q[3];
rz(-1.67098) q[3];
sx q[3];
rz(-1.8377931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8835939) q[0];
sx q[0];
rz(-0.086253919) q[0];
sx q[0];
rz(-1.8121383) q[0];
rz(-2.6138002) q[1];
sx q[1];
rz(-1.3161491) q[1];
sx q[1];
rz(-2.7425308) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51080158) q[0];
sx q[0];
rz(-0.36040655) q[0];
sx q[0];
rz(-1.7921757) q[0];
x q[1];
rz(-1.9827323) q[2];
sx q[2];
rz(-0.88232458) q[2];
sx q[2];
rz(0.43718034) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0001155) q[1];
sx q[1];
rz(-0.85073227) q[1];
sx q[1];
rz(1.1517555) q[1];
rz(-1.7979223) q[3];
sx q[3];
rz(-2.4339614) q[3];
sx q[3];
rz(-2.9766102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.993678) q[2];
sx q[2];
rz(-0.75239158) q[2];
sx q[2];
rz(-1.0047151) q[2];
rz(-0.29804722) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(1.0684439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3598969) q[0];
sx q[0];
rz(-1.5157461) q[0];
sx q[0];
rz(-1.9500649) q[0];
rz(0.53939348) q[1];
sx q[1];
rz(-2.0627779) q[1];
sx q[1];
rz(0.76621145) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82220413) q[0];
sx q[0];
rz(-2.4766716) q[0];
sx q[0];
rz(-2.1524236) q[0];
x q[1];
rz(2.2257588) q[2];
sx q[2];
rz(-2.1127491) q[2];
sx q[2];
rz(2.4905027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6886544) q[1];
sx q[1];
rz(-2.3804745) q[1];
sx q[1];
rz(-1.903694) q[1];
rz(-pi) q[2];
rz(-0.096682354) q[3];
sx q[3];
rz(-1.9396693) q[3];
sx q[3];
rz(2.3547805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3300276) q[2];
sx q[2];
rz(-1.6946239) q[2];
sx q[2];
rz(-2.2583029) q[2];
rz(-0.53089321) q[3];
sx q[3];
rz(-0.58404946) q[3];
sx q[3];
rz(-2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74060488) q[0];
sx q[0];
rz(-0.29380909) q[0];
sx q[0];
rz(0.64939943) q[0];
rz(-0.3736639) q[1];
sx q[1];
rz(-2.3039736) q[1];
sx q[1];
rz(1.1423473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.861688) q[0];
sx q[0];
rz(-1.7293469) q[0];
sx q[0];
rz(-1.9025361) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5552005) q[2];
sx q[2];
rz(-1.7806539) q[2];
sx q[2];
rz(3.1210528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1159276) q[1];
sx q[1];
rz(-0.84995228) q[1];
sx q[1];
rz(-1.9864818) q[1];
x q[2];
rz(1.1116696) q[3];
sx q[3];
rz(-1.5304426) q[3];
sx q[3];
rz(-0.99276517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1335527) q[2];
sx q[2];
rz(-2.4456761) q[2];
sx q[2];
rz(-0.35888654) q[2];
rz(-2.533203) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(-0.86148328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1818039) q[0];
sx q[0];
rz(-2.602808) q[0];
sx q[0];
rz(-0.96610075) q[0];
rz(0.44202647) q[1];
sx q[1];
rz(-1.288488) q[1];
sx q[1];
rz(2.7332222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.934487) q[0];
sx q[0];
rz(-1.7033332) q[0];
sx q[0];
rz(0.026959628) q[0];
rz(-pi) q[1];
rz(2.9938575) q[2];
sx q[2];
rz(-1.6134173) q[2];
sx q[2];
rz(-0.10382596) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4736259) q[1];
sx q[1];
rz(-1.7136116) q[1];
sx q[1];
rz(-2.1599139) q[1];
rz(-0.73159742) q[3];
sx q[3];
rz(-1.7642191) q[3];
sx q[3];
rz(0.0091656589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24892204) q[2];
sx q[2];
rz(-2.9327124) q[2];
sx q[2];
rz(1.199031) q[2];
rz(-1.8438953) q[3];
sx q[3];
rz(-1.0782995) q[3];
sx q[3];
rz(-0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93160981) q[0];
sx q[0];
rz(-1.4298121) q[0];
sx q[0];
rz(-2.0344875) q[0];
rz(-2.9529086) q[1];
sx q[1];
rz(-2.76177) q[1];
sx q[1];
rz(1.3245378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8134091) q[0];
sx q[0];
rz(-1.8556917) q[0];
sx q[0];
rz(-2.6064742) q[0];
rz(0.021090551) q[2];
sx q[2];
rz(-2.167835) q[2];
sx q[2];
rz(1.326438) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75211084) q[1];
sx q[1];
rz(-0.5403203) q[1];
sx q[1];
rz(-1.575927) q[1];
rz(-pi) q[2];
rz(-1.0774433) q[3];
sx q[3];
rz(-2.2098118) q[3];
sx q[3];
rz(2.8324156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.79335) q[2];
sx q[2];
rz(-2.5620354) q[2];
sx q[2];
rz(0.14774518) q[2];
rz(-0.82792264) q[3];
sx q[3];
rz(-1.4180276) q[3];
sx q[3];
rz(2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70656908) q[0];
sx q[0];
rz(-0.35025418) q[0];
sx q[0];
rz(-1.6178004) q[0];
rz(-1.2500457) q[1];
sx q[1];
rz(-0.83386546) q[1];
sx q[1];
rz(-2.7203383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4604302) q[0];
sx q[0];
rz(-1.958485) q[0];
sx q[0];
rz(1.9676542) q[0];
rz(-1.658861) q[2];
sx q[2];
rz(-0.95092623) q[2];
sx q[2];
rz(1.6234835) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4688465) q[1];
sx q[1];
rz(-1.4839659) q[1];
sx q[1];
rz(-1.561291) q[1];
rz(2.6700745) q[3];
sx q[3];
rz(-0.49905159) q[3];
sx q[3];
rz(1.7744482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0201575) q[2];
sx q[2];
rz(-2.406106) q[2];
sx q[2];
rz(-2.6911733) q[2];
rz(-1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(-2.569017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12238518) q[0];
sx q[0];
rz(-1.5337802) q[0];
sx q[0];
rz(1.5564729) q[0];
rz(2.4284594) q[1];
sx q[1];
rz(-1.5819555) q[1];
sx q[1];
rz(-3.1176007) q[1];
rz(1.0456741) q[2];
sx q[2];
rz(-1.9782421) q[2];
sx q[2];
rz(-1.1876456) q[2];
rz(-1.6433257) q[3];
sx q[3];
rz(-1.4909496) q[3];
sx q[3];
rz(-0.95076233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
