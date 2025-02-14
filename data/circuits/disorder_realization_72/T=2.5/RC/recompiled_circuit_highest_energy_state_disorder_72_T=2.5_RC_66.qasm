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
rz(0.82062757) q[0];
sx q[0];
rz(-0.66250116) q[0];
sx q[0];
rz(2.878046) q[0];
rz(1.1072371) q[1];
sx q[1];
rz(4.4237408) q[1];
sx q[1];
rz(10.10034) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1206664) q[0];
sx q[0];
rz(-1.5772737) q[0];
sx q[0];
rz(1.9064039) q[0];
rz(-0.58568875) q[2];
sx q[2];
rz(-2.0721451) q[2];
sx q[2];
rz(-0.96282178) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61900422) q[1];
sx q[1];
rz(-1.1824885) q[1];
sx q[1];
rz(0.16242811) q[1];
rz(-2.1249902) q[3];
sx q[3];
rz(-0.57376353) q[3];
sx q[3];
rz(-0.91780182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2519309) q[2];
sx q[2];
rz(-1.2222922) q[2];
sx q[2];
rz(2.7737854) q[2];
rz(2.8273072) q[3];
sx q[3];
rz(-2.8511484) q[3];
sx q[3];
rz(0.17816003) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3322068) q[0];
sx q[0];
rz(-0.089501373) q[0];
sx q[0];
rz(-1.3595164) q[0];
rz(1.8571732) q[1];
sx q[1];
rz(-1.9764683) q[1];
sx q[1];
rz(3.0899835) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415539) q[0];
sx q[0];
rz(-1.7518028) q[0];
sx q[0];
rz(1.8405528) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9412291) q[2];
sx q[2];
rz(-1.7631754) q[2];
sx q[2];
rz(-2.7967601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1882798) q[1];
sx q[1];
rz(-1.9447535) q[1];
sx q[1];
rz(0.00061051173) q[1];
x q[2];
rz(-1.1647928) q[3];
sx q[3];
rz(-2.0116802) q[3];
sx q[3];
rz(-2.1390381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.08903683) q[2];
sx q[2];
rz(-0.26036981) q[2];
sx q[2];
rz(0.54711771) q[2];
rz(-1.268092) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(-0.33856302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38076213) q[0];
sx q[0];
rz(-1.0365423) q[0];
sx q[0];
rz(1.3408252) q[0];
rz(2.2876168) q[1];
sx q[1];
rz(-1.8546591) q[1];
sx q[1];
rz(-0.30240789) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.069347) q[0];
sx q[0];
rz(-2.4503142) q[0];
sx q[0];
rz(-2.012558) q[0];
rz(-2.688681) q[2];
sx q[2];
rz(-2.6476417) q[2];
sx q[2];
rz(3.0978705) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21061347) q[1];
sx q[1];
rz(-0.96081464) q[1];
sx q[1];
rz(0.90175866) q[1];
x q[2];
rz(-0.27389483) q[3];
sx q[3];
rz(-1.0948101) q[3];
sx q[3];
rz(2.7416503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2053908) q[2];
sx q[2];
rz(-1.8855636) q[2];
sx q[2];
rz(0.54226792) q[2];
rz(2.1451456) q[3];
sx q[3];
rz(-0.22765972) q[3];
sx q[3];
rz(-3.0136133) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19313136) q[0];
sx q[0];
rz(-1.4840935) q[0];
sx q[0];
rz(-1.2633854) q[0];
rz(0.43426934) q[1];
sx q[1];
rz(-2.3691005) q[1];
sx q[1];
rz(1.6901406) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5853441) q[0];
sx q[0];
rz(-0.55216575) q[0];
sx q[0];
rz(0.57741965) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0201597) q[2];
sx q[2];
rz(-1.0673293) q[2];
sx q[2];
rz(-1.2502647) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5401104) q[1];
sx q[1];
rz(-2.9992691) q[1];
sx q[1];
rz(-2.2156634) q[1];
rz(0.55813467) q[3];
sx q[3];
rz(-1.0252748) q[3];
sx q[3];
rz(-1.6203875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0666535) q[2];
sx q[2];
rz(-0.94330072) q[2];
sx q[2];
rz(-2.2554876) q[2];
rz(3.1352299) q[3];
sx q[3];
rz(-1.8013835) q[3];
sx q[3];
rz(-2.0224915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5804382) q[0];
sx q[0];
rz(-2.6851324) q[0];
sx q[0];
rz(1.891267) q[0];
rz(2.9087032) q[1];
sx q[1];
rz(-2.3857375) q[1];
sx q[1];
rz(-0.09045352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1250153) q[0];
sx q[0];
rz(-0.27864405) q[0];
sx q[0];
rz(-2.0296627) q[0];
rz(-pi) q[1];
rz(-1.301342) q[2];
sx q[2];
rz(-1.1134992) q[2];
sx q[2];
rz(-2.2895165) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4793746) q[1];
sx q[1];
rz(-1.5081128) q[1];
sx q[1];
rz(-2.4880243) q[1];
rz(-0.085898475) q[3];
sx q[3];
rz(-1.2441564) q[3];
sx q[3];
rz(2.896691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.49454409) q[2];
sx q[2];
rz(-2.9477951) q[2];
sx q[2];
rz(-1.2079283) q[2];
rz(0.9211933) q[3];
sx q[3];
rz(-2.2974206) q[3];
sx q[3];
rz(-2.6827961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10814609) q[0];
sx q[0];
rz(-1.7195846) q[0];
sx q[0];
rz(-2.6500927) q[0];
rz(-0.72832251) q[1];
sx q[1];
rz(-2.0305384) q[1];
sx q[1];
rz(2.947015) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74516836) q[0];
sx q[0];
rz(-2.130742) q[0];
sx q[0];
rz(0.26225175) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75677432) q[2];
sx q[2];
rz(-0.96454731) q[2];
sx q[2];
rz(-0.9204677) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7214365) q[1];
sx q[1];
rz(-2.2682574) q[1];
sx q[1];
rz(-2.2583794) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65690144) q[3];
sx q[3];
rz(-1.3852775) q[3];
sx q[3];
rz(2.3168062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3538094) q[2];
sx q[2];
rz(-2.1259191) q[2];
sx q[2];
rz(3.1397528) q[2];
rz(-1.6483824) q[3];
sx q[3];
rz(-2.4108678) q[3];
sx q[3];
rz(0.28436896) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581185) q[0];
sx q[0];
rz(-0.26171568) q[0];
sx q[0];
rz(-2.0765685) q[0];
rz(0.26271543) q[1];
sx q[1];
rz(-2.5698667) q[1];
sx q[1];
rz(-0.46331847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36998591) q[0];
sx q[0];
rz(-0.60877234) q[0];
sx q[0];
rz(-1.2941525) q[0];
rz(-2.5538374) q[2];
sx q[2];
rz(-0.81771933) q[2];
sx q[2];
rz(0.14434926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.66760705) q[1];
sx q[1];
rz(-1.4346448) q[1];
sx q[1];
rz(-2.6546247) q[1];
rz(1.7780532) q[3];
sx q[3];
rz(-0.1175783) q[3];
sx q[3];
rz(-2.1685947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.00039438417) q[2];
sx q[2];
rz(-1.980915) q[2];
sx q[2];
rz(-2.311643) q[2];
rz(-2.0937894) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(-0.1786264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038427) q[0];
sx q[0];
rz(-2.1512845) q[0];
sx q[0];
rz(-2.3867699) q[0];
rz(-2.7524475) q[1];
sx q[1];
rz(-1.6837589) q[1];
sx q[1];
rz(2.7438502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6104113) q[0];
sx q[0];
rz(-1.8540148) q[0];
sx q[0];
rz(0.2679268) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5903274) q[2];
sx q[2];
rz(-1.3756075) q[2];
sx q[2];
rz(-0.27182394) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6036247) q[1];
sx q[1];
rz(-1.84314) q[1];
sx q[1];
rz(-0.26428278) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92614545) q[3];
sx q[3];
rz(-2.20813) q[3];
sx q[3];
rz(-1.4755352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8336746) q[2];
sx q[2];
rz(-0.29611823) q[2];
sx q[2];
rz(-2.8274242) q[2];
rz(1.2416154) q[3];
sx q[3];
rz(-1.5528468) q[3];
sx q[3];
rz(1.2392905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.817713) q[0];
sx q[0];
rz(-0.68809026) q[0];
sx q[0];
rz(2.896198) q[0];
rz(1.9112401) q[1];
sx q[1];
rz(-0.10086682) q[1];
sx q[1];
rz(2.6255677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300872) q[0];
sx q[0];
rz(-1.1908414) q[0];
sx q[0];
rz(-1.7854693) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86271) q[2];
sx q[2];
rz(-2.1893775) q[2];
sx q[2];
rz(-1.0363621) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5616331) q[1];
sx q[1];
rz(-1.206534) q[1];
sx q[1];
rz(-1.8749401) q[1];
rz(-1.5088169) q[3];
sx q[3];
rz(-1.7070642) q[3];
sx q[3];
rz(-2.800022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55769354) q[2];
sx q[2];
rz(-1.2217174) q[2];
sx q[2];
rz(-1.0830797) q[2];
rz(0.22100581) q[3];
sx q[3];
rz(-0.14328863) q[3];
sx q[3];
rz(2.3451282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021127064) q[0];
sx q[0];
rz(-0.080862232) q[0];
sx q[0];
rz(3.0057111) q[0];
rz(1.7120301) q[1];
sx q[1];
rz(-1.1212965) q[1];
sx q[1];
rz(0.28724614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7571952) q[0];
sx q[0];
rz(-1.383019) q[0];
sx q[0];
rz(-1.3223935) q[0];
rz(-pi) q[1];
rz(0.084566074) q[2];
sx q[2];
rz(-1.8075602) q[2];
sx q[2];
rz(-2.6344476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8441287) q[1];
sx q[1];
rz(-1.7904568) q[1];
sx q[1];
rz(2.541887) q[1];
rz(-pi) q[2];
rz(-2.1549018) q[3];
sx q[3];
rz(-1.3249448) q[3];
sx q[3];
rz(0.81013269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2222432) q[2];
sx q[2];
rz(-1.2390077) q[2];
sx q[2];
rz(-2.3075721) q[2];
rz(-0.30759865) q[3];
sx q[3];
rz(-2.7643876) q[3];
sx q[3];
rz(2.8789177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1766227) q[0];
sx q[0];
rz(-1.8360092) q[0];
sx q[0];
rz(-1.0259884) q[0];
rz(-0.38656825) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(-2.8564524) q[2];
sx q[2];
rz(-2.7013473) q[2];
sx q[2];
rz(-0.57821748) q[2];
rz(-0.14561226) q[3];
sx q[3];
rz(-1.3619193) q[3];
sx q[3];
rz(-3.0633434) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
