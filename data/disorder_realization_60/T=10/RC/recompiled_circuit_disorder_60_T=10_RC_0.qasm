OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.94937593) q[0];
sx q[0];
rz(5.2360143) q[0];
sx q[0];
rz(9.4935023) q[0];
rz(1.7460495) q[1];
sx q[1];
rz(-1.6092602) q[1];
sx q[1];
rz(-1.2083763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3291408) q[0];
sx q[0];
rz(-1.666823) q[0];
sx q[0];
rz(-3.0949233) q[0];
x q[1];
rz(-1.5829093) q[2];
sx q[2];
rz(-1.6965869) q[2];
sx q[2];
rz(2.8319401) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0304058) q[1];
sx q[1];
rz(-0.62796794) q[1];
sx q[1];
rz(-2.9524132) q[1];
rz(-0.73322202) q[3];
sx q[3];
rz(-0.65922046) q[3];
sx q[3];
rz(2.1472907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8157114) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.6750083) q[2];
rz(0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(0.74716032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.117347) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(-2.970447) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(0.29719621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3575465) q[0];
sx q[0];
rz(-1.6539126) q[0];
sx q[0];
rz(3.0931285) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4291123) q[2];
sx q[2];
rz(-1.2006294) q[2];
sx q[2];
rz(2.9287101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5479996) q[1];
sx q[1];
rz(-1.9108692) q[1];
sx q[1];
rz(2.1151357) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76241775) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(-3.1114515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.979636) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(-2.5276108) q[2];
rz(-2.2654514) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0959594) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(2.8339548) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(-0.83980733) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1727985) q[0];
sx q[0];
rz(-1.286924) q[0];
sx q[0];
rz(2.7127405) q[0];
rz(0.33939056) q[2];
sx q[2];
rz(-1.4815785) q[2];
sx q[2];
rz(2.3238376) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3629344) q[1];
sx q[1];
rz(-2.0211453) q[1];
sx q[1];
rz(1.910701) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93785502) q[3];
sx q[3];
rz(-1.4301392) q[3];
sx q[3];
rz(1.3682287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86747375) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(1.3179368) q[2];
rz(1.2157724) q[3];
sx q[3];
rz(-0.35651818) q[3];
sx q[3];
rz(1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777305) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(-2.5207649) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(-0.96558085) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25528204) q[0];
sx q[0];
rz(-0.83034407) q[0];
sx q[0];
rz(-1.1925973) q[0];
rz(-pi) q[1];
rz(3.1104286) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(-0.63809168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5406148) q[1];
sx q[1];
rz(-1.317273) q[1];
sx q[1];
rz(-0.032932245) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2772917) q[3];
sx q[3];
rz(-2.9923277) q[3];
sx q[3];
rz(2.5316558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3199557) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(-1.2840282) q[0];
rz(2.8517826) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(-1.0481542) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6402733) q[0];
sx q[0];
rz(-0.87448705) q[0];
sx q[0];
rz(0.93080824) q[0];
x q[1];
rz(-2.3088147) q[2];
sx q[2];
rz(-1.1156429) q[2];
sx q[2];
rz(2.8684794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54742868) q[1];
sx q[1];
rz(-2.8058726) q[1];
sx q[1];
rz(0.28433849) q[1];
rz(-0.9976451) q[3];
sx q[3];
rz(-0.76612681) q[3];
sx q[3];
rz(-2.0875974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(2.2407545) q[2];
rz(2.0488996) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(1.2341011) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(2.545488) q[0];
rz(-1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.8966819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0743474) q[0];
sx q[0];
rz(-1.0676358) q[0];
sx q[0];
rz(0.36672451) q[0];
rz(-pi) q[1];
rz(1.2539358) q[2];
sx q[2];
rz(-2.0358026) q[2];
sx q[2];
rz(0.76380619) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5925496) q[1];
sx q[1];
rz(-2.2144496) q[1];
sx q[1];
rz(-2.2085269) q[1];
x q[2];
rz(-2.8935029) q[3];
sx q[3];
rz(-1.3458985) q[3];
sx q[3];
rz(3.0450862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.231455) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(2.9166252) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(-0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769619) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(0.27994573) q[0];
rz(-1.4631368) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(-2.8889012) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1278909) q[0];
sx q[0];
rz(-1.2780006) q[0];
sx q[0];
rz(-0.049833628) q[0];
rz(-3.1227038) q[2];
sx q[2];
rz(-0.9364555) q[2];
sx q[2];
rz(2.6332476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6934769) q[1];
sx q[1];
rz(-0.87065334) q[1];
sx q[1];
rz(-0.6853939) q[1];
x q[2];
rz(-2.3971862) q[3];
sx q[3];
rz(-0.91074569) q[3];
sx q[3];
rz(0.50124121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8213886) q[2];
sx q[2];
rz(-0.89367047) q[2];
sx q[2];
rz(1.5318711) q[2];
rz(-1.948471) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(-1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(3.0367592) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(2.8727818) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(0.2789467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9439518) q[0];
sx q[0];
rz(-2.1924207) q[0];
sx q[0];
rz(2.5191514) q[0];
x q[1];
rz(-2.9697044) q[2];
sx q[2];
rz(-1.0749146) q[2];
sx q[2];
rz(-1.9672729) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5178691) q[1];
sx q[1];
rz(-1.1964799) q[1];
sx q[1];
rz(-2.3831297) q[1];
rz(-1.7989743) q[3];
sx q[3];
rz(-0.5842714) q[3];
sx q[3];
rz(-0.97464558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(1.0507978) q[3];
sx q[3];
rz(-1.2069586) q[3];
sx q[3];
rz(-1.1999493) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779697) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(1.8883702) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(-1.995283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86664591) q[0];
sx q[0];
rz(-1.5005611) q[0];
sx q[0];
rz(-2.3932891) q[0];
rz(-pi) q[1];
rz(-0.65593221) q[2];
sx q[2];
rz(-2.1428875) q[2];
sx q[2];
rz(-0.42665542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4790926) q[1];
sx q[1];
rz(-0.88516419) q[1];
sx q[1];
rz(-0.094866026) q[1];
rz(-pi) q[2];
rz(0.54543145) q[3];
sx q[3];
rz(-0.71119961) q[3];
sx q[3];
rz(-1.4382854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15381947) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(-2.231797) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290454) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(3.1273499) q[0];
rz(2.2968538) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(-1.3815809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68122278) q[0];
sx q[0];
rz(-1.5813706) q[0];
sx q[0];
rz(-1.4282385) q[0];
rz(-pi) q[1];
rz(0.99190418) q[2];
sx q[2];
rz(-1.5169946) q[2];
sx q[2];
rz(-2.3827041) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33243079) q[1];
sx q[1];
rz(-0.5756439) q[1];
sx q[1];
rz(-1.8368506) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53346177) q[3];
sx q[3];
rz(-1.7603612) q[3];
sx q[3];
rz(1.5272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0766729) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(0.89407095) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0513231) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(1.8016215) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(-2.4253035) q[2];
sx q[2];
rz(-1.204797) q[2];
sx q[2];
rz(2.7439678) q[2];
rz(-2.1469231) q[3];
sx q[3];
rz(-2.4352286) q[3];
sx q[3];
rz(2.2027204) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
