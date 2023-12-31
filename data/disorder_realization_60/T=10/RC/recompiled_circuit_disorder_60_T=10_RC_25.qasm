OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1922167) q[0];
sx q[0];
rz(-2.0944216) q[0];
sx q[0];
rz(3.0728683) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(1.2083763) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24613334) q[0];
sx q[0];
rz(-1.6172505) q[0];
sx q[0];
rz(1.4746656) q[0];
rz(-3.0460998) q[2];
sx q[2];
rz(-3.0152233) q[2];
sx q[2];
rz(-0.40590826) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5282643) q[1];
sx q[1];
rz(-1.4600888) q[1];
sx q[1];
rz(0.61943357) q[1];
rz(2.049202) q[3];
sx q[3];
rz(-1.0983101) q[3];
sx q[3];
rz(-0.14379584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8157114) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.6750083) q[2];
rz(0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(-2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0242457) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(-1.1741937) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(-0.29719621) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8288119) q[0];
sx q[0];
rz(-0.096185616) q[0];
sx q[0];
rz(-2.0975153) q[0];
x q[1];
rz(-0.71248033) q[2];
sx q[2];
rz(-1.2006294) q[2];
sx q[2];
rz(-0.21288255) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5479996) q[1];
sx q[1];
rz(-1.9108692) q[1];
sx q[1];
rz(1.026457) q[1];
x q[2];
rz(-1.9678715) q[3];
sx q[3];
rz(-2.2928772) q[3];
sx q[3];
rz(1.811036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(-0.87614122) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(2.8339548) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(2.3017853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6162286) q[0];
sx q[0];
rz(-1.9814241) q[0];
sx q[0];
rz(1.2603659) q[0];
rz(-pi) q[1];
rz(-1.4762127) q[2];
sx q[2];
rz(-1.9087831) q[2];
sx q[2];
rz(0.72159492) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0452022) q[1];
sx q[1];
rz(-0.55711105) q[1];
sx q[1];
rz(-2.5379009) q[1];
x q[2];
rz(-2.2037376) q[3];
sx q[3];
rz(-1.7114534) q[3];
sx q[3];
rz(1.7733639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(1.8236558) q[2];
rz(1.9258202) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(-1.6320451) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76386219) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(2.6960301) q[0];
rz(-0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(2.1760118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27698101) q[0];
sx q[0];
rz(-0.81482139) q[0];
sx q[0];
rz(-2.7576202) q[0];
rz(-pi) q[1];
rz(-0.84237174) q[2];
sx q[2];
rz(-1.5475376) q[2];
sx q[2];
rz(-2.1881441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5406148) q[1];
sx q[1];
rz(-1.317273) q[1];
sx q[1];
rz(3.1086604) q[1];
x q[2];
rz(-1.4569034) q[3];
sx q[3];
rz(-1.4741065) q[3];
sx q[3];
rz(0.25988042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.821637) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(2.2122673) q[2];
rz(2.4980513) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(-1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(1.8575645) q[0];
rz(2.8517826) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(2.0934385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832344) q[0];
sx q[0];
rz(-0.90792197) q[0];
sx q[0];
rz(-2.5213084) q[0];
rz(-pi) q[1];
rz(-0.94190188) q[2];
sx q[2];
rz(-2.2976029) q[2];
sx q[2];
rz(0.84743365) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2939799) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(-1.4732248) q[1];
x q[2];
rz(0.48091472) q[3];
sx q[3];
rz(-2.1926583) q[3];
sx q[3];
rz(-2.8180168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-0.90083814) q[2];
rz(-1.0926931) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(1.9074915) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95935217) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(2.545488) q[0];
rz(-1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4620074) q[0];
sx q[0];
rz(-1.8903268) q[0];
sx q[0];
rz(-1.0380448) q[0];
rz(2.5858324) q[2];
sx q[2];
rz(-2.5854977) q[2];
sx q[2];
rz(-0.13242002) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60331261) q[1];
sx q[1];
rz(-1.0744175) q[1];
sx q[1];
rz(0.75116317) q[1];
rz(-pi) q[2];
rz(0.75002589) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(-2.1961574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.231455) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(-2.9166252) q[2];
rz(-0.088430017) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769619) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(0.27994573) q[0];
rz(1.4631368) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(-0.25269145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7130816) q[0];
sx q[0];
rz(-1.6185074) q[0];
sx q[0];
rz(-1.2776572) q[0];
rz(-pi) q[1];
rz(-2.2052223) q[2];
sx q[2];
rz(-1.5860103) q[2];
sx q[2];
rz(-1.051256) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.44811571) q[1];
sx q[1];
rz(-0.87065334) q[1];
sx q[1];
rz(2.4561988) q[1];
x q[2];
rz(-2.3971862) q[3];
sx q[3];
rz(-0.91074569) q[3];
sx q[3];
rz(0.50124121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8213886) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(1.6097216) q[2];
rz(1.1931217) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(-1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0367592) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(2.8727818) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(2.862646) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3725961) q[0];
sx q[0];
rz(-1.0770174) q[0];
sx q[0];
rz(2.2934224) q[0];
x q[1];
rz(-0.17188822) q[2];
sx q[2];
rz(-1.0749146) q[2];
sx q[2];
rz(-1.1743197) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3153509) q[1];
sx q[1];
rz(-0.8289753) q[1];
sx q[1];
rz(0.51893236) q[1];
rz(-pi) q[2];
rz(2.9931077) q[3];
sx q[3];
rz(-1.0035702) q[3];
sx q[3];
rz(0.70311577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(-1.5489244) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46362296) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(-1.2532225) q[0];
rz(0.62250096) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(1.1463096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2749467) q[0];
sx q[0];
rz(-1.6410315) q[0];
sx q[0];
rz(-0.74830351) q[0];
rz(-pi) q[1];
rz(-2.253138) q[2];
sx q[2];
rz(-1.0323712) q[2];
sx q[2];
rz(0.74935645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.031530023) q[1];
sx q[1];
rz(-1.6441802) q[1];
sx q[1];
rz(-2.258638) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5067234) q[3];
sx q[3];
rz(-1.2253309) q[3];
sx q[3];
rz(-2.8431901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9877732) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(-2.1257607) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(2.773496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(3.1273499) q[0];
rz(-0.8447389) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(1.7600118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.325557) q[0];
sx q[0];
rz(-2.9986458) q[0];
sx q[0];
rz(-1.4965034) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99190418) q[2];
sx q[2];
rz(-1.5169946) q[2];
sx q[2];
rz(-0.75888854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8091619) q[1];
sx q[1];
rz(-0.5756439) q[1];
sx q[1];
rz(-1.8368506) q[1];
rz(-1.3515477) q[3];
sx q[3];
rz(-1.0478813) q[3];
sx q[3];
rz(0.067283665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0766729) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(2.1515965) q[2];
rz(0.89407095) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(-0.22542424) q[3];
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
rz(-pi/2) q[3];
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
rz(-1.0513231) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(1.8016215) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(-1.1006533) q[2];
sx q[2];
rz(-2.2307776) q[2];
sx q[2];
rz(-2.270436) q[2];
rz(2.1918478) q[3];
sx q[3];
rz(-1.2093778) q[3];
sx q[3];
rz(-2.9686684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
