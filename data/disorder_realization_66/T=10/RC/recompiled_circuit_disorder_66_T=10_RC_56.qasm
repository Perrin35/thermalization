OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(1.6605641) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(3.0552157) q[1];
sx q[1];
rz(9.4019158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829553) q[0];
sx q[0];
rz(-1.4541172) q[0];
sx q[0];
rz(-1.4215683) q[0];
rz(-pi) q[1];
rz(2.0055662) q[2];
sx q[2];
rz(-2.044894) q[2];
sx q[2];
rz(0.20027645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7658246) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(-1.9181812) q[1];
x q[2];
rz(-2.1894737) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(-1.4116956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(3.0549333) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21762411) q[0];
sx q[0];
rz(-2.1930165) q[0];
sx q[0];
rz(-1.5909965) q[0];
rz(-pi) q[1];
rz(2.8367963) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(3.0457029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3955235) q[1];
sx q[1];
rz(-1.6568526) q[1];
sx q[1];
rz(1.6834016) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3195023) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(3.086123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(1.6764199) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6572606) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(1.3844301) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4577643) q[0];
sx q[0];
rz(-1.0965075) q[0];
sx q[0];
rz(-1.477924) q[0];
rz(1.0648107) q[2];
sx q[2];
rz(-1.2227321) q[2];
sx q[2];
rz(0.51007523) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82921925) q[1];
sx q[1];
rz(-0.81554283) q[1];
sx q[1];
rz(-1.6014598) q[1];
x q[2];
rz(-2.3178187) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(-2.2658474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8973792) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(-1.3440075) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(0.20733325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509388) q[0];
sx q[0];
rz(-1.5293984) q[0];
sx q[0];
rz(0.008965094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7447128) q[2];
sx q[2];
rz(-1.5408278) q[2];
sx q[2];
rz(-2.7001691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9878848) q[1];
sx q[1];
rz(-2.8578836) q[1];
sx q[1];
rz(-0.99516408) q[1];
rz(2.6299423) q[3];
sx q[3];
rz(-1.5582065) q[3];
sx q[3];
rz(0.97012855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31071445) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(-0.80777848) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(-2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(-0.89548683) q[0];
rz(-0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(-0.53363824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69913188) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(-1.2752227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76782121) q[2];
sx q[2];
rz(-1.4066833) q[2];
sx q[2];
rz(2.5196911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.071760885) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(-1.7204309) q[1];
x q[2];
rz(0.20547262) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(-1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(-2.1369381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6001119) q[0];
sx q[0];
rz(-1.2280557) q[0];
sx q[0];
rz(-0.39370763) q[0];
rz(-pi) q[1];
rz(-2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(-2.4900988) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5914454) q[1];
sx q[1];
rz(-2.1131556) q[1];
sx q[1];
rz(-1.28246) q[1];
rz(-pi) q[2];
rz(-1.1569571) q[3];
sx q[3];
rz(-1.5598179) q[3];
sx q[3];
rz(1.699284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7375609) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(-1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(0.67214322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37049609) q[0];
sx q[0];
rz(-1.5633977) q[0];
sx q[0];
rz(-1.7367944) q[0];
rz(-pi) q[1];
rz(2.3357046) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(-2.9758331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0969095) q[1];
sx q[1];
rz(-1.9367366) q[1];
sx q[1];
rz(2.6198322) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70127212) q[3];
sx q[3];
rz(-1.4605195) q[3];
sx q[3];
rz(2.9993204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44052112) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-2.3349082) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8950302) q[0];
sx q[0];
rz(-2.6916305) q[0];
sx q[0];
rz(-1.297784) q[0];
x q[1];
rz(-2.9358747) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(-0.92668698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17864171) q[1];
sx q[1];
rz(-0.9045524) q[1];
sx q[1];
rz(1.2013749) q[1];
rz(-pi) q[2];
rz(1.673647) q[3];
sx q[3];
rz(-0.30234435) q[3];
sx q[3];
rz(-0.24149382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(0.65752423) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-0.34067571) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9690902) q[0];
sx q[0];
rz(-1.889125) q[0];
sx q[0];
rz(0.83842917) q[0];
x q[1];
rz(2.9195243) q[2];
sx q[2];
rz(-1.0169612) q[2];
sx q[2];
rz(2.1332707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.454168) q[1];
sx q[1];
rz(-1.7523645) q[1];
sx q[1];
rz(-2.1329692) q[1];
rz(0.19823234) q[3];
sx q[3];
rz(-1.4054338) q[3];
sx q[3];
rz(-0.73393047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72934812) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-0.27858946) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(0.07671193) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4519276) q[0];
sx q[0];
rz(-1.8849012) q[0];
sx q[0];
rz(3.0118224) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3086583) q[2];
sx q[2];
rz(-0.58912504) q[2];
sx q[2];
rz(2.4208456) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1288651) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(-0.058458316) q[1];
rz(2.6807908) q[3];
sx q[3];
rz(-2.9446903) q[3];
sx q[3];
rz(2.868696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(3.07807) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615622) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(0.2858328) q[2];
sx q[2];
rz(-1.1163859) q[2];
sx q[2];
rz(-3.1039539) q[2];
rz(0.42701957) q[3];
sx q[3];
rz(-1.9058766) q[3];
sx q[3];
rz(1.8591892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
