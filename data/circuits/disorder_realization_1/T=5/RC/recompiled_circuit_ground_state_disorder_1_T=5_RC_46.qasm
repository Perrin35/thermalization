OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27788568) q[0];
sx q[0];
rz(-2.7122893) q[0];
sx q[0];
rz(2.8306146) q[0];
rz(1.7929945) q[1];
sx q[1];
rz(-1.0962948) q[1];
sx q[1];
rz(-2.5675093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35168655) q[0];
sx q[0];
rz(-1.154019) q[0];
sx q[0];
rz(0.4051286) q[0];
rz(-0.034717807) q[2];
sx q[2];
rz(-0.98811921) q[2];
sx q[2];
rz(2.5947844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8694386) q[1];
sx q[1];
rz(-1.660214) q[1];
sx q[1];
rz(1.2050259) q[1];
rz(3.0776377) q[3];
sx q[3];
rz(-1.1950781) q[3];
sx q[3];
rz(2.2967754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2744039) q[2];
sx q[2];
rz(-1.2588661) q[2];
sx q[2];
rz(-1.3192419) q[2];
rz(0.073171767) q[3];
sx q[3];
rz(-1.8354445) q[3];
sx q[3];
rz(0.81899548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0480334) q[0];
sx q[0];
rz(-1.9608542) q[0];
sx q[0];
rz(2.4617885) q[0];
rz(-2.3381084) q[1];
sx q[1];
rz(-1.3619224) q[1];
sx q[1];
rz(-2.5018073) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8713952) q[0];
sx q[0];
rz(-1.2551487) q[0];
sx q[0];
rz(-0.90243602) q[0];
rz(-2.4117208) q[2];
sx q[2];
rz(-1.6449071) q[2];
sx q[2];
rz(0.48219901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3181139) q[1];
sx q[1];
rz(-0.9645071) q[1];
sx q[1];
rz(2.4068474) q[1];
x q[2];
rz(-1.3469668) q[3];
sx q[3];
rz(-1.8719684) q[3];
sx q[3];
rz(2.6909242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5378319) q[2];
sx q[2];
rz(-2.64309) q[2];
sx q[2];
rz(2.4375088) q[2];
rz(-3.0294561) q[3];
sx q[3];
rz(-1.4487368) q[3];
sx q[3];
rz(-2.1191547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10244399) q[0];
sx q[0];
rz(-1.1071438) q[0];
sx q[0];
rz(0.33945864) q[0];
rz(-2.0231694) q[1];
sx q[1];
rz(-0.92670852) q[1];
sx q[1];
rz(-0.035645398) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4887071) q[0];
sx q[0];
rz(-1.6566601) q[0];
sx q[0];
rz(-0.20767943) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72398846) q[2];
sx q[2];
rz(-0.44000726) q[2];
sx q[2];
rz(-2.2429446) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0782203) q[1];
sx q[1];
rz(-2.2148892) q[1];
sx q[1];
rz(-1.3011342) q[1];
x q[2];
rz(-2.6958353) q[3];
sx q[3];
rz(-0.66346079) q[3];
sx q[3];
rz(-1.9395246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3493335) q[2];
sx q[2];
rz(-2.4534241) q[2];
sx q[2];
rz(2.2815857) q[2];
rz(2.3490014) q[3];
sx q[3];
rz(-2.9383797) q[3];
sx q[3];
rz(2.121076) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.129313) q[0];
sx q[0];
rz(-1.7565933) q[0];
sx q[0];
rz(2.5087575) q[0];
rz(1.9889779) q[1];
sx q[1];
rz(-0.8747789) q[1];
sx q[1];
rz(-1.4702183) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7790079) q[0];
sx q[0];
rz(-0.92155313) q[0];
sx q[0];
rz(0.016223793) q[0];
x q[1];
rz(-0.15369065) q[2];
sx q[2];
rz(-2.3665303) q[2];
sx q[2];
rz(0.77010051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51293514) q[1];
sx q[1];
rz(-2.3875065) q[1];
sx q[1];
rz(-1.9505461) q[1];
rz(-pi) q[2];
rz(0.65514647) q[3];
sx q[3];
rz(-2.0787424) q[3];
sx q[3];
rz(-0.20482132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4319438) q[2];
sx q[2];
rz(-0.96082965) q[2];
sx q[2];
rz(0.43080899) q[2];
rz(2.4591947) q[3];
sx q[3];
rz(-1.6513446) q[3];
sx q[3];
rz(-2.1916913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7948941) q[0];
sx q[0];
rz(-1.8759202) q[0];
sx q[0];
rz(-1.1899813) q[0];
rz(2.8311912) q[1];
sx q[1];
rz(-2.0605395) q[1];
sx q[1];
rz(0.50500542) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9146945) q[0];
sx q[0];
rz(-2.3932308) q[0];
sx q[0];
rz(2.2415461) q[0];
rz(-pi) q[1];
rz(1.9857384) q[2];
sx q[2];
rz(-1.6908852) q[2];
sx q[2];
rz(-1.1121554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47520275) q[1];
sx q[1];
rz(-2.0670927) q[1];
sx q[1];
rz(1.6700891) q[1];
x q[2];
rz(-1.4851863) q[3];
sx q[3];
rz(-2.2946649) q[3];
sx q[3];
rz(-1.6023265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7348822) q[2];
sx q[2];
rz(-2.4411185) q[2];
sx q[2];
rz(-0.90275466) q[2];
rz(-1.0550176) q[3];
sx q[3];
rz(-2.2746634) q[3];
sx q[3];
rz(-2.6796851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91586739) q[0];
sx q[0];
rz(-0.71284717) q[0];
sx q[0];
rz(-2.0676887) q[0];
rz(1.7621) q[1];
sx q[1];
rz(-0.72934377) q[1];
sx q[1];
rz(-3.0440547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9363234) q[0];
sx q[0];
rz(-2.5714834) q[0];
sx q[0];
rz(0.3400712) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9768841) q[2];
sx q[2];
rz(-0.40702074) q[2];
sx q[2];
rz(-0.87813745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8733858) q[1];
sx q[1];
rz(-1.6318351) q[1];
sx q[1];
rz(0.048442099) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0076526) q[3];
sx q[3];
rz(-1.1112257) q[3];
sx q[3];
rz(-0.28467049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11233106) q[2];
sx q[2];
rz(-1.6882201) q[2];
sx q[2];
rz(2.7304999) q[2];
rz(1.7279651) q[3];
sx q[3];
rz(-2.3634383) q[3];
sx q[3];
rz(1.5467862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77808648) q[0];
sx q[0];
rz(-3.111105) q[0];
sx q[0];
rz(-1.4084858) q[0];
rz(-2.1259437) q[1];
sx q[1];
rz(-1.8135095) q[1];
sx q[1];
rz(1.6633063) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41970872) q[0];
sx q[0];
rz(-2.1745178) q[0];
sx q[0];
rz(-1.9323856) q[0];
x q[1];
rz(0.64817102) q[2];
sx q[2];
rz(-1.0888466) q[2];
sx q[2];
rz(0.75343695) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46131733) q[1];
sx q[1];
rz(-1.6467386) q[1];
sx q[1];
rz(1.126775) q[1];
x q[2];
rz(2.0709345) q[3];
sx q[3];
rz(-1.6864417) q[3];
sx q[3];
rz(-1.0341737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0279072) q[2];
sx q[2];
rz(-0.62782851) q[2];
sx q[2];
rz(0.17024635) q[2];
rz(1.2804735) q[3];
sx q[3];
rz(-1.6672641) q[3];
sx q[3];
rz(1.1580275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4729507) q[0];
sx q[0];
rz(-0.96422115) q[0];
sx q[0];
rz(-0.52841312) q[0];
rz(0.93327418) q[1];
sx q[1];
rz(-2.2163138) q[1];
sx q[1];
rz(-2.5859213) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816102) q[0];
sx q[0];
rz(-0.9564119) q[0];
sx q[0];
rz(-0.23991628) q[0];
rz(-pi) q[1];
rz(-3.0425983) q[2];
sx q[2];
rz(-0.79455355) q[2];
sx q[2];
rz(1.0459448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.56951504) q[1];
sx q[1];
rz(-1.6028499) q[1];
sx q[1];
rz(0.7116913) q[1];
x q[2];
rz(2.7714588) q[3];
sx q[3];
rz(-2.4151093) q[3];
sx q[3];
rz(-2.671482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8133424) q[2];
sx q[2];
rz(-0.80208653) q[2];
sx q[2];
rz(0.30725202) q[2];
rz(-0.79636374) q[3];
sx q[3];
rz(-0.73085228) q[3];
sx q[3];
rz(-3.0439607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55427134) q[0];
sx q[0];
rz(-2.0926496) q[0];
sx q[0];
rz(1.3294504) q[0];
rz(-0.21996552) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(1.3185917) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585478) q[0];
sx q[0];
rz(-1.3038583) q[0];
sx q[0];
rz(0.22375317) q[0];
rz(-pi) q[1];
rz(-1.4542142) q[2];
sx q[2];
rz(-2.0676842) q[2];
sx q[2];
rz(-1.0036482) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.66178759) q[1];
sx q[1];
rz(-2.0347682) q[1];
sx q[1];
rz(-0.3825835) q[1];
x q[2];
rz(-0.72033003) q[3];
sx q[3];
rz(-1.2033278) q[3];
sx q[3];
rz(-3.0331963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4465176) q[2];
sx q[2];
rz(-1.7138285) q[2];
sx q[2];
rz(1.8565149) q[2];
rz(-0.88589969) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(1.6133962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1000243) q[0];
sx q[0];
rz(-2.4055241) q[0];
sx q[0];
rz(2.8247483) q[0];
rz(2.7007804) q[1];
sx q[1];
rz(-0.79439729) q[1];
sx q[1];
rz(0.37030927) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99890868) q[0];
sx q[0];
rz(-2.7591191) q[0];
sx q[0];
rz(-1.6539025) q[0];
rz(-2.2787408) q[2];
sx q[2];
rz(-0.86332317) q[2];
sx q[2];
rz(1.6528636) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.036344254) q[1];
sx q[1];
rz(-2.7351028) q[1];
sx q[1];
rz(-2.2986733) q[1];
rz(0.050906128) q[3];
sx q[3];
rz(-1.6744782) q[3];
sx q[3];
rz(-0.10242505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4050498) q[2];
sx q[2];
rz(-2.6506347) q[2];
sx q[2];
rz(-1.6323818) q[2];
rz(-2.3368808) q[3];
sx q[3];
rz(-1.5038306) q[3];
sx q[3];
rz(3.0338083) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0263154) q[0];
sx q[0];
rz(-1.4125217) q[0];
sx q[0];
rz(-1.9052292) q[0];
rz(2.9112877) q[1];
sx q[1];
rz(-2.8254012) q[1];
sx q[1];
rz(0.91513035) q[1];
rz(0.66441734) q[2];
sx q[2];
rz(-1.0894486) q[2];
sx q[2];
rz(1.1846381) q[2];
rz(-2.0143916) q[3];
sx q[3];
rz(-0.66772912) q[3];
sx q[3];
rz(-2.276788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
