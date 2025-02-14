OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11256448) q[0];
sx q[0];
rz(-1.636314) q[0];
sx q[0];
rz(-2.3583052) q[0];
rz(-2.0593491) q[1];
sx q[1];
rz(4.6286197) q[1];
sx q[1];
rz(11.774465) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5675674) q[0];
sx q[0];
rz(-1.8301395) q[0];
sx q[0];
rz(-1.7631084) q[0];
rz(-pi) q[1];
rz(-1.7246805) q[2];
sx q[2];
rz(-1.6072011) q[2];
sx q[2];
rz(-2.3652181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4152849) q[1];
sx q[1];
rz(-1.0276762) q[1];
sx q[1];
rz(1.5341137) q[1];
x q[2];
rz(0.57609419) q[3];
sx q[3];
rz(-0.8818501) q[3];
sx q[3];
rz(0.16932936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5504494) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(3.0421416) q[2];
rz(-0.24672306) q[3];
sx q[3];
rz(-1.5244502) q[3];
sx q[3];
rz(-0.39768404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0196911) q[0];
sx q[0];
rz(-2.1653403) q[0];
sx q[0];
rz(-0.9285399) q[0];
rz(-1.6909201) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(-2.4931152) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6715393) q[0];
sx q[0];
rz(-2.3190365) q[0];
sx q[0];
rz(2.3683192) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5669539) q[2];
sx q[2];
rz(-1.840072) q[2];
sx q[2];
rz(0.93967162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14997031) q[1];
sx q[1];
rz(-2.9987965) q[1];
sx q[1];
rz(-3.0317329) q[1];
x q[2];
rz(-3.0233211) q[3];
sx q[3];
rz(-1.857548) q[3];
sx q[3];
rz(2.3210382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4497946) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(-0.85255426) q[2];
rz(-2.2954588) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(2.1107296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344675) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(-2.0563828) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(1.6403713) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1000993) q[0];
sx q[0];
rz(-1.3736808) q[0];
sx q[0];
rz(1.6105349) q[0];
rz(0.88884647) q[2];
sx q[2];
rz(-0.29817981) q[2];
sx q[2];
rz(-2.2928638) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5021421) q[1];
sx q[1];
rz(-1.0002563) q[1];
sx q[1];
rz(-3.023663) q[1];
x q[2];
rz(-0.15764938) q[3];
sx q[3];
rz(-1.356989) q[3];
sx q[3];
rz(0.42663867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7368855) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(2.286818) q[2];
rz(-2.2331623) q[3];
sx q[3];
rz(-1.5769438) q[3];
sx q[3];
rz(-1.1165718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.214355) q[0];
sx q[0];
rz(-0.4011918) q[0];
sx q[0];
rz(1.9450564) q[0];
rz(-1.475097) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(1.1845142) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9036983) q[0];
sx q[0];
rz(-1.8082976) q[0];
sx q[0];
rz(-0.71126513) q[0];
rz(-pi) q[1];
rz(-1.3802211) q[2];
sx q[2];
rz(-0.82697502) q[2];
sx q[2];
rz(2.8883052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6946751) q[1];
sx q[1];
rz(-1.0694587) q[1];
sx q[1];
rz(0.25872177) q[1];
x q[2];
rz(-0.27359815) q[3];
sx q[3];
rz(-1.7471385) q[3];
sx q[3];
rz(0.30418744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2148332) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(-1.2614177) q[2];
rz(0.43241209) q[3];
sx q[3];
rz(-2.0093446) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686907) q[0];
sx q[0];
rz(-1.9111159) q[0];
sx q[0];
rz(3.1121837) q[0];
rz(-0.75621653) q[1];
sx q[1];
rz(-2.5543946) q[1];
sx q[1];
rz(-2.3888033) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2761951) q[0];
sx q[0];
rz(-1.5711693) q[0];
sx q[0];
rz(-1.5764109) q[0];
x q[1];
rz(1.3026313) q[2];
sx q[2];
rz(-1.6577621) q[2];
sx q[2];
rz(3.027107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1323586) q[1];
sx q[1];
rz(-1.0028603) q[1];
sx q[1];
rz(-1.6949937) q[1];
rz(-2.4642508) q[3];
sx q[3];
rz(-0.65653446) q[3];
sx q[3];
rz(1.7650106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0011562) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(-0.39503869) q[2];
rz(-0.47518528) q[3];
sx q[3];
rz(-1.7893712) q[3];
sx q[3];
rz(-2.7648259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.7659371) q[0];
sx q[0];
rz(-2.4212615) q[0];
sx q[0];
rz(0.96250594) q[0];
rz(2.6267701) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(0.33822507) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048934919) q[0];
sx q[0];
rz(-0.83507628) q[0];
sx q[0];
rz(-0.2661163) q[0];
x q[1];
rz(1.8184616) q[2];
sx q[2];
rz(-1.3405352) q[2];
sx q[2];
rz(2.3375653) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6392614) q[1];
sx q[1];
rz(-3.1384472) q[1];
sx q[1];
rz(1.5622713) q[1];
x q[2];
rz(-0.068840222) q[3];
sx q[3];
rz(-1.1496115) q[3];
sx q[3];
rz(-2.0023458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6414791) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(-2.5872453) q[2];
rz(2.2902299) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(0.62801492) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8165269) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(0.75041962) q[0];
rz(0.57506192) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(0.65779984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3391425) q[0];
sx q[0];
rz(-0.93451148) q[0];
sx q[0];
rz(2.0872981) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5267685) q[2];
sx q[2];
rz(-2.2775473) q[2];
sx q[2];
rz(-2.6296774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21322528) q[1];
sx q[1];
rz(-1.2847273) q[1];
sx q[1];
rz(-0.017315344) q[1];
rz(-pi) q[2];
rz(-2.0714893) q[3];
sx q[3];
rz(-0.4288396) q[3];
sx q[3];
rz(-2.2718112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.647992) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(-1.6592525) q[2];
rz(3.0209387) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(-0.83546662) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7928612) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(0.29801512) q[0];
rz(1.0385849) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(0.15377741) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6385495) q[0];
sx q[0];
rz(-0.26493236) q[0];
sx q[0];
rz(0.28980906) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91570274) q[2];
sx q[2];
rz(-0.88353744) q[2];
sx q[2];
rz(2.5926431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.41339918) q[1];
sx q[1];
rz(-1.6604042) q[1];
sx q[1];
rz(2.7941568) q[1];
x q[2];
rz(1.1438391) q[3];
sx q[3];
rz(-1.4597963) q[3];
sx q[3];
rz(0.75253651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0457354) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.9528961) q[2];
rz(1.3304322) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(0.73330283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7981912) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(1.4991624) q[0];
rz(-2.5921953) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(2.8909491) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17535398) q[0];
sx q[0];
rz(-0.73161517) q[0];
sx q[0];
rz(-2.0425955) q[0];
rz(-pi) q[1];
x q[1];
rz(1.556384) q[2];
sx q[2];
rz(-1.8317128) q[2];
sx q[2];
rz(-3.122807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.033015164) q[1];
sx q[1];
rz(-0.77627722) q[1];
sx q[1];
rz(-3.0423711) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6969243) q[3];
sx q[3];
rz(-1.3775004) q[3];
sx q[3];
rz(-2.2910558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9833019) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(1.0164725) q[2];
rz(-0.086183444) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(-2.3763954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8330399) q[0];
sx q[0];
rz(-2.2873531) q[0];
sx q[0];
rz(0.41900751) q[0];
rz(2.6047756) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(-2.4057665) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2013071) q[0];
sx q[0];
rz(-2.2129411) q[0];
sx q[0];
rz(1.8117732) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0334098) q[2];
sx q[2];
rz(-0.99796406) q[2];
sx q[2];
rz(2.9705641) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6169301) q[1];
sx q[1];
rz(-1.7866275) q[1];
sx q[1];
rz(-0.642435) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9189776) q[3];
sx q[3];
rz(-1.1135808) q[3];
sx q[3];
rz(1.3570076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1824823) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(2.518892) q[2];
rz(-0.73838082) q[3];
sx q[3];
rz(-1.1850971) q[3];
sx q[3];
rz(-0.27899376) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644792) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(2.4979757) q[1];
sx q[1];
rz(-1.2871965) q[1];
sx q[1];
rz(1.0866477) q[1];
rz(1.0254597) q[2];
sx q[2];
rz(-0.42517107) q[2];
sx q[2];
rz(-0.33381768) q[2];
rz(-0.52538659) q[3];
sx q[3];
rz(-0.27670866) q[3];
sx q[3];
rz(1.8075734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
