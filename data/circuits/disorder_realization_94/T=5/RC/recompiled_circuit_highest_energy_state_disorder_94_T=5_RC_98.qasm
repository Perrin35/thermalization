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
rz(0.87297451) q[0];
sx q[0];
rz(2.8277446) q[0];
sx q[0];
rz(10.477973) q[0];
rz(-1.516951) q[1];
sx q[1];
rz(-2.8487974) q[1];
sx q[1];
rz(2.204978) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40799403) q[0];
sx q[0];
rz(-2.6331365) q[0];
sx q[0];
rz(-1.5336799) q[0];
rz(-0.83990015) q[2];
sx q[2];
rz(-1.1625625) q[2];
sx q[2];
rz(-2.3917835) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9642311) q[1];
sx q[1];
rz(-1.291591) q[1];
sx q[1];
rz(-2.5482168) q[1];
rz(-pi) q[2];
rz(-1.7630834) q[3];
sx q[3];
rz(-1.1858127) q[3];
sx q[3];
rz(-0.92038233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8006353) q[2];
sx q[2];
rz(-2.1850047) q[2];
sx q[2];
rz(2.0802278) q[2];
rz(2.495885) q[3];
sx q[3];
rz(-0.3568477) q[3];
sx q[3];
rz(2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2467932) q[0];
sx q[0];
rz(-2.6428887) q[0];
sx q[0];
rz(-2.6820768) q[0];
rz(-1.3872321) q[1];
sx q[1];
rz(-2.6315755) q[1];
sx q[1];
rz(1.0646819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014039847) q[0];
sx q[0];
rz(-1.8888374) q[0];
sx q[0];
rz(-1.7377338) q[0];
rz(-1.5534723) q[2];
sx q[2];
rz(-1.8380062) q[2];
sx q[2];
rz(1.653275) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.83512703) q[1];
sx q[1];
rz(-1.4773932) q[1];
sx q[1];
rz(0.3208092) q[1];
rz(-pi) q[2];
rz(-2.2467733) q[3];
sx q[3];
rz(-2.7193391) q[3];
sx q[3];
rz(0.32988068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.77949828) q[2];
sx q[2];
rz(-0.67216122) q[2];
sx q[2];
rz(-0.54074311) q[2];
rz(2.8703459) q[3];
sx q[3];
rz(-1.2416779) q[3];
sx q[3];
rz(2.0103683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3830477) q[0];
sx q[0];
rz(-0.31065148) q[0];
sx q[0];
rz(-2.7591822) q[0];
rz(-0.27579871) q[1];
sx q[1];
rz(-2.661992) q[1];
sx q[1];
rz(-0.84024215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0226757) q[0];
sx q[0];
rz(-1.2166326) q[0];
sx q[0];
rz(-2.7849331) q[0];
x q[1];
rz(-2.1539168) q[2];
sx q[2];
rz(-0.35735574) q[2];
sx q[2];
rz(0.59229702) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0663755) q[1];
sx q[1];
rz(-1.0569968) q[1];
sx q[1];
rz(-2.2877076) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48583416) q[3];
sx q[3];
rz(-1.5155025) q[3];
sx q[3];
rz(-3.0818617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8648839) q[2];
sx q[2];
rz(-1.0397006) q[2];
sx q[2];
rz(0.060982171) q[2];
rz(1.3649155) q[3];
sx q[3];
rz(-2.958332) q[3];
sx q[3];
rz(2.0257559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404496) q[0];
sx q[0];
rz(-1.018486) q[0];
sx q[0];
rz(-2.2818991) q[0];
rz(2.1172093) q[1];
sx q[1];
rz(-0.59737098) q[1];
sx q[1];
rz(0.76622564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47536248) q[0];
sx q[0];
rz(-1.6138617) q[0];
sx q[0];
rz(-1.5353139) q[0];
rz(1.0417074) q[2];
sx q[2];
rz(-1.7211564) q[2];
sx q[2];
rz(2.9632115) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1658881) q[1];
sx q[1];
rz(-2.5448826) q[1];
sx q[1];
rz(2.6590682) q[1];
rz(-3.0697508) q[3];
sx q[3];
rz(-0.86962442) q[3];
sx q[3];
rz(-1.5506684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1102981) q[2];
sx q[2];
rz(-0.86250192) q[2];
sx q[2];
rz(3.1150225) q[2];
rz(-0.26818141) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(0.69494438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7705695) q[0];
sx q[0];
rz(-0.71582782) q[0];
sx q[0];
rz(-0.41719607) q[0];
rz(-0.5873276) q[1];
sx q[1];
rz(-2.6081577) q[1];
sx q[1];
rz(2.3951098) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6339109) q[0];
sx q[0];
rz(-0.48134781) q[0];
sx q[0];
rz(-2.4018628) q[0];
x q[1];
rz(0.76144371) q[2];
sx q[2];
rz(-0.9599182) q[2];
sx q[2];
rz(2.2100984) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7938817) q[1];
sx q[1];
rz(-1.0679967) q[1];
sx q[1];
rz(-2.0050843) q[1];
rz(-1.0899441) q[3];
sx q[3];
rz(-2.0084511) q[3];
sx q[3];
rz(-0.57695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.31009659) q[2];
sx q[2];
rz(-0.3336755) q[2];
sx q[2];
rz(2.2200072) q[2];
rz(-2.0690252) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(0.5630365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801341) q[0];
sx q[0];
rz(-0.38668329) q[0];
sx q[0];
rz(2.4899546) q[0];
rz(2.1561275) q[1];
sx q[1];
rz(-2.5990504) q[1];
sx q[1];
rz(-0.14269565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7008377) q[0];
sx q[0];
rz(-1.6336226) q[0];
sx q[0];
rz(1.2962893) q[0];
x q[1];
rz(1.5161523) q[2];
sx q[2];
rz(-1.3846591) q[2];
sx q[2];
rz(0.7344377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.922077) q[1];
sx q[1];
rz(-1.0990881) q[1];
sx q[1];
rz(-1.707915) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0246932) q[3];
sx q[3];
rz(-1.4479234) q[3];
sx q[3];
rz(-1.529913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5720713) q[2];
sx q[2];
rz(-0.60939747) q[2];
sx q[2];
rz(-0.93192464) q[2];
rz(2.7052687) q[3];
sx q[3];
rz(-2.6659129) q[3];
sx q[3];
rz(-1.1599734) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6465004) q[0];
sx q[0];
rz(-0.91630542) q[0];
sx q[0];
rz(1.6012993) q[0];
rz(1.0013927) q[1];
sx q[1];
rz(-1.5903558) q[1];
sx q[1];
rz(-0.95759773) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0315899) q[0];
sx q[0];
rz(-2.6837807) q[0];
sx q[0];
rz(-0.97719595) q[0];
rz(-pi) q[1];
rz(0.25156542) q[2];
sx q[2];
rz(-1.6645164) q[2];
sx q[2];
rz(2.9523015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55064978) q[1];
sx q[1];
rz(-0.47524449) q[1];
sx q[1];
rz(-1.5859423) q[1];
rz(-1.6041338) q[3];
sx q[3];
rz(-2.1289724) q[3];
sx q[3];
rz(-2.6128925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5074978) q[2];
sx q[2];
rz(-1.8562506) q[2];
sx q[2];
rz(1.1320587) q[2];
rz(3.0105528) q[3];
sx q[3];
rz(-1.9877501) q[3];
sx q[3];
rz(-1.2426144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89192724) q[0];
sx q[0];
rz(-0.51814336) q[0];
sx q[0];
rz(-3.0666572) q[0];
rz(-2.9451008) q[1];
sx q[1];
rz(-0.34759977) q[1];
sx q[1];
rz(-0.67952716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.644548) q[0];
sx q[0];
rz(-1.0168075) q[0];
sx q[0];
rz(-1.333167) q[0];
x q[1];
rz(1.1668233) q[2];
sx q[2];
rz(-1.5979287) q[2];
sx q[2];
rz(-0.89744842) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7767002) q[1];
sx q[1];
rz(-2.6235995) q[1];
sx q[1];
rz(-2.1692976) q[1];
rz(0.51107589) q[3];
sx q[3];
rz(-1.5835488) q[3];
sx q[3];
rz(-1.2619357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3139265) q[2];
sx q[2];
rz(-0.9845261) q[2];
sx q[2];
rz(1.0796245) q[2];
rz(0.45352724) q[3];
sx q[3];
rz(-0.87117666) q[3];
sx q[3];
rz(-2.1760904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51665783) q[0];
sx q[0];
rz(-2.9624532) q[0];
sx q[0];
rz(3.0058885) q[0];
rz(-0.16595674) q[1];
sx q[1];
rz(-0.9981007) q[1];
sx q[1];
rz(-0.96806324) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0128202) q[0];
sx q[0];
rz(-0.51542437) q[0];
sx q[0];
rz(1.5933871) q[0];
rz(-2.1130457) q[2];
sx q[2];
rz(-1.442277) q[2];
sx q[2];
rz(-0.1947386) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9688473) q[1];
sx q[1];
rz(-2.0683534) q[1];
sx q[1];
rz(1.8109591) q[1];
x q[2];
rz(-1.8289205) q[3];
sx q[3];
rz(-1.2069697) q[3];
sx q[3];
rz(2.8175333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3193937) q[2];
sx q[2];
rz(-0.92331702) q[2];
sx q[2];
rz(-0.53259069) q[2];
rz(-2.3432664) q[3];
sx q[3];
rz(-0.39755487) q[3];
sx q[3];
rz(3.047191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3548729) q[0];
sx q[0];
rz(-0.69364554) q[0];
sx q[0];
rz(-0.68674809) q[0];
rz(1.2314388) q[1];
sx q[1];
rz(-1.310692) q[1];
sx q[1];
rz(0.062006921) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3162295) q[0];
sx q[0];
rz(-2.4333409) q[0];
sx q[0];
rz(-3.0407719) q[0];
rz(-pi) q[1];
rz(2.7786534) q[2];
sx q[2];
rz(-0.20359765) q[2];
sx q[2];
rz(2.3516977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7534291) q[1];
sx q[1];
rz(-1.7723188) q[1];
sx q[1];
rz(2.9934197) q[1];
x q[2];
rz(2.8436321) q[3];
sx q[3];
rz(-1.5567257) q[3];
sx q[3];
rz(1.5009106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9318781) q[2];
sx q[2];
rz(-2.176602) q[2];
sx q[2];
rz(-2.9689201) q[2];
rz(2.4058345) q[3];
sx q[3];
rz(-2.5685205) q[3];
sx q[3];
rz(0.47634038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.944735) q[0];
sx q[0];
rz(-1.4516964) q[0];
sx q[0];
rz(-1.4781937) q[0];
rz(0.860515) q[1];
sx q[1];
rz(-1.9445226) q[1];
sx q[1];
rz(-1.3938211) q[1];
rz(-1.6780268) q[2];
sx q[2];
rz(-1.6926395) q[2];
sx q[2];
rz(-2.0933685) q[2];
rz(-2.2205856) q[3];
sx q[3];
rz(-2.3719127) q[3];
sx q[3];
rz(1.0815892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
