OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(-1.4927827) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(-0.69505039) q[1];
sx q[1];
rz(1.0746497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5045835) q[0];
sx q[0];
rz(-1.552581) q[0];
sx q[0];
rz(1.5760742) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4104112) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(0.10057848) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1638012) q[1];
sx q[1];
rz(-1.4492387) q[1];
sx q[1];
rz(-1.6519283) q[1];
rz(-2.0243458) q[3];
sx q[3];
rz(-2.0587066) q[3];
sx q[3];
rz(-0.32353668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.6281698) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56869498) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.5011903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26186865) q[0];
sx q[0];
rz(-0.12517087) q[0];
sx q[0];
rz(0.20010389) q[0];
rz(-pi) q[1];
rz(-1.00374) q[2];
sx q[2];
rz(-0.54076414) q[2];
sx q[2];
rz(-0.31976779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3850296) q[1];
sx q[1];
rz(-0.90365138) q[1];
sx q[1];
rz(1.7900311) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5562015) q[3];
sx q[3];
rz(-2.7505034) q[3];
sx q[3];
rz(0.99816834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8555277) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(1.1668011) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353772) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(1.0789543) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.5037781) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8684981) q[0];
sx q[0];
rz(-1.6582186) q[0];
sx q[0];
rz(-3.0762061) q[0];
rz(-2.2321646) q[2];
sx q[2];
rz(-2.4436908) q[2];
sx q[2];
rz(2.1991889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40239247) q[1];
sx q[1];
rz(-1.270073) q[1];
sx q[1];
rz(-1.6154358) q[1];
x q[2];
rz(2.2272439) q[3];
sx q[3];
rz(-1.4667061) q[3];
sx q[3];
rz(-0.57746938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7109795) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(2.3977996) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1733615) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(-0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7950492) q[0];
sx q[0];
rz(-1.8115037) q[0];
sx q[0];
rz(1.4562777) q[0];
rz(-pi) q[1];
rz(-0.97781424) q[2];
sx q[2];
rz(-1.2056418) q[2];
sx q[2];
rz(-2.1821373) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60050868) q[1];
sx q[1];
rz(-0.87227548) q[1];
sx q[1];
rz(-2.7545616) q[1];
x q[2];
rz(-1.9150919) q[3];
sx q[3];
rz(-0.85840423) q[3];
sx q[3];
rz(0.32574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(2.4728298) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(-2.0892129) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(0.88422424) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28750557) q[0];
sx q[0];
rz(-0.26212087) q[0];
sx q[0];
rz(0.8192807) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82489478) q[2];
sx q[2];
rz(-0.56606228) q[2];
sx q[2];
rz(-2.538946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21806949) q[1];
sx q[1];
rz(-1.4265718) q[1];
sx q[1];
rz(-1.6378535) q[1];
rz(-pi) q[2];
rz(2.7719284) q[3];
sx q[3];
rz(-2.0737994) q[3];
sx q[3];
rz(-0.30077416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9050682) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(-2.0489676) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(-2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(-2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1554402) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(-0.040779671) q[0];
x q[1];
rz(0.71408748) q[2];
sx q[2];
rz(-1.519515) q[2];
sx q[2];
rz(2.4209792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5507257) q[1];
sx q[1];
rz(-2.496965) q[1];
sx q[1];
rz(-2.4025737) q[1];
rz(-pi) q[2];
rz(1.045741) q[3];
sx q[3];
rz(-1.9805555) q[3];
sx q[3];
rz(2.4214782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(0.058241025) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(1.4490022) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.70599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5933701) q[0];
sx q[0];
rz(-2.3832294) q[0];
sx q[0];
rz(2.3286657) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5510467) q[2];
sx q[2];
rz(-1.5524459) q[2];
sx q[2];
rz(1.8576647) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.095852921) q[1];
sx q[1];
rz(-2.6624655) q[1];
sx q[1];
rz(-0.87289255) q[1];
rz(-1.0649879) q[3];
sx q[3];
rz(-1.9516264) q[3];
sx q[3];
rz(1.8246458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-2.984916) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85132861) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(3.0086349) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(-1.7763604) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93850219) q[0];
sx q[0];
rz(-1.7068958) q[0];
sx q[0];
rz(0.076119856) q[0];
x q[1];
rz(-1.4012785) q[2];
sx q[2];
rz(-3.0119544) q[2];
sx q[2];
rz(0.6946176) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3268765) q[1];
sx q[1];
rz(-1.5183581) q[1];
sx q[1];
rz(-2.609842) q[1];
rz(-pi) q[2];
rz(-2.7216464) q[3];
sx q[3];
rz(-1.4010324) q[3];
sx q[3];
rz(1.1743634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.075295538) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(-2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(-0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(2.2176567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16238427) q[0];
sx q[0];
rz(-0.93063336) q[0];
sx q[0];
rz(1.0248653) q[0];
rz(-pi) q[1];
rz(0.48859476) q[2];
sx q[2];
rz(-1.5813507) q[2];
sx q[2];
rz(0.34305629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82211557) q[1];
sx q[1];
rz(-1.176782) q[1];
sx q[1];
rz(-1.0072717) q[1];
x q[2];
rz(-0.71420788) q[3];
sx q[3];
rz(-1.076477) q[3];
sx q[3];
rz(0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.3941992) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(0.53846255) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(-1.3884397) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-0.038169233) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725806) q[0];
sx q[0];
rz(-0.60927143) q[0];
sx q[0];
rz(-1.6174181) q[0];
rz(-pi) q[1];
x q[1];
rz(1.795904) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(3.0181146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8074179) q[1];
sx q[1];
rz(-1.5290698) q[1];
sx q[1];
rz(1.9445022) q[1];
rz(-pi) q[2];
rz(0.60572259) q[3];
sx q[3];
rz(-2.1401792) q[3];
sx q[3];
rz(-0.13525087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.9434628) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643628) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(-2.6565068) q[2];
sx q[2];
rz(-3.017364) q[2];
sx q[2];
rz(-0.45807522) q[2];
rz(-1.2002767) q[3];
sx q[3];
rz(-0.73642086) q[3];
sx q[3];
rz(0.45263276) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
