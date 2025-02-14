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
rz(3.570896) q[0];
sx q[0];
rz(12.255393) q[0];
rz(1.7929945) q[1];
sx q[1];
rz(-1.0962948) q[1];
sx q[1];
rz(-2.5675093) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0943756) q[0];
sx q[0];
rz(-1.939491) q[0];
sx q[0];
rz(-1.1218907) q[0];
rz(-pi) q[1];
rz(1.5181731) q[2];
sx q[2];
rz(-0.5835909) q[2];
sx q[2];
rz(0.4837732) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2644538) q[1];
sx q[1];
rz(-1.2065556) q[1];
sx q[1];
rz(0.095714494) q[1];
rz(1.4101683) q[3];
sx q[3];
rz(-0.38086748) q[3];
sx q[3];
rz(-2.4695652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8671888) q[2];
sx q[2];
rz(-1.2588661) q[2];
sx q[2];
rz(1.8223507) q[2];
rz(-0.073171767) q[3];
sx q[3];
rz(-1.8354445) q[3];
sx q[3];
rz(-0.81899548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.0935593) q[0];
sx q[0];
rz(-1.9608542) q[0];
sx q[0];
rz(-0.67980415) q[0];
rz(0.80348429) q[1];
sx q[1];
rz(-1.3619224) q[1];
sx q[1];
rz(-2.5018073) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6006192) q[0];
sx q[0];
rz(-2.2007211) q[0];
sx q[0];
rz(0.39430228) q[0];
rz(-pi) q[1];
rz(1.670094) q[2];
sx q[2];
rz(-2.298215) q[2];
sx q[2];
rz(1.0224487) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8315994) q[1];
sx q[1];
rz(-2.2265456) q[1];
sx q[1];
rz(2.3393231) q[1];
rz(-0.30840318) q[3];
sx q[3];
rz(-1.7843912) q[3];
sx q[3];
rz(-1.9540389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5378319) q[2];
sx q[2];
rz(-2.64309) q[2];
sx q[2];
rz(-0.70408386) q[2];
rz(0.1121366) q[3];
sx q[3];
rz(-1.6928558) q[3];
sx q[3];
rz(-1.0224379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10244399) q[0];
sx q[0];
rz(-2.0344489) q[0];
sx q[0];
rz(-2.802134) q[0];
rz(-1.1184232) q[1];
sx q[1];
rz(-2.2148841) q[1];
sx q[1];
rz(3.1059473) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6729926) q[0];
sx q[0];
rz(-0.22449271) q[0];
sx q[0];
rz(-2.7461282) q[0];
rz(-pi) q[1];
rz(-0.33907922) q[2];
sx q[2];
rz(-1.8568175) q[2];
sx q[2];
rz(-1.3468483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0782203) q[1];
sx q[1];
rz(-0.92670346) q[1];
sx q[1];
rz(1.8404585) q[1];
x q[2];
rz(-0.44575739) q[3];
sx q[3];
rz(-0.66346079) q[3];
sx q[3];
rz(1.9395246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3493335) q[2];
sx q[2];
rz(-0.68816853) q[2];
sx q[2];
rz(2.2815857) q[2];
rz(0.79259121) q[3];
sx q[3];
rz(-0.20321295) q[3];
sx q[3];
rz(-1.0205166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.129313) q[0];
sx q[0];
rz(-1.3849994) q[0];
sx q[0];
rz(-2.5087575) q[0];
rz(1.9889779) q[1];
sx q[1];
rz(-2.2668138) q[1];
sx q[1];
rz(-1.6713743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8058384) q[0];
sx q[0];
rz(-2.4921761) q[0];
sx q[0];
rz(1.592167) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4219513) q[2];
sx q[2];
rz(-2.334377) q[2];
sx q[2];
rz(-2.5850353) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51293514) q[1];
sx q[1];
rz(-0.75408616) q[1];
sx q[1];
rz(1.1910466) q[1];
x q[2];
rz(-2.4864462) q[3];
sx q[3];
rz(-1.0628502) q[3];
sx q[3];
rz(0.20482132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4319438) q[2];
sx q[2];
rz(-2.180763) q[2];
sx q[2];
rz(2.7107837) q[2];
rz(0.68239799) q[3];
sx q[3];
rz(-1.4902481) q[3];
sx q[3];
rz(-2.1916913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7948941) q[0];
sx q[0];
rz(-1.8759202) q[0];
sx q[0];
rz(-1.9516113) q[0];
rz(-2.8311912) q[1];
sx q[1];
rz(-2.0605395) q[1];
sx q[1];
rz(-0.50500542) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22689817) q[0];
sx q[0];
rz(-0.74836181) q[0];
sx q[0];
rz(-2.2415461) q[0];
x q[1];
rz(1.9857384) q[2];
sx q[2];
rz(-1.6908852) q[2];
sx q[2];
rz(2.0294373) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4601535) q[1];
sx q[1];
rz(-2.6362754) q[1];
sx q[1];
rz(-0.1810592) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4851863) q[3];
sx q[3];
rz(-2.2946649) q[3];
sx q[3];
rz(-1.5392661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7348822) q[2];
sx q[2];
rz(-0.70047417) q[2];
sx q[2];
rz(-0.90275466) q[2];
rz(-2.086575) q[3];
sx q[3];
rz(-2.2746634) q[3];
sx q[3];
rz(2.6796851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91586739) q[0];
sx q[0];
rz(-0.71284717) q[0];
sx q[0];
rz(1.073904) q[0];
rz(-1.7621) q[1];
sx q[1];
rz(-0.72934377) q[1];
sx q[1];
rz(-0.097537907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65501761) q[0];
sx q[0];
rz(-1.3897822) q[0];
sx q[0];
rz(0.54365309) q[0];
x q[1];
rz(-2.9729207) q[2];
sx q[2];
rz(-1.198581) q[2];
sx q[2];
rz(-1.8255359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3055468) q[1];
sx q[1];
rz(-1.5224445) q[1];
sx q[1];
rz(1.6319066) q[1];
rz(-0.52953665) q[3];
sx q[3];
rz(-2.0696928) q[3];
sx q[3];
rz(-1.0130817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11233106) q[2];
sx q[2];
rz(-1.4533726) q[2];
sx q[2];
rz(-0.41109273) q[2];
rz(1.4136275) q[3];
sx q[3];
rz(-0.77815431) q[3];
sx q[3];
rz(1.5467862) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3635062) q[0];
sx q[0];
rz(-0.030487617) q[0];
sx q[0];
rz(1.4084858) q[0];
rz(2.1259437) q[1];
sx q[1];
rz(-1.3280832) q[1];
sx q[1];
rz(1.6633063) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7218839) q[0];
sx q[0];
rz(-2.1745178) q[0];
sx q[0];
rz(1.9323856) q[0];
rz(-0.71395771) q[2];
sx q[2];
rz(-0.78642008) q[2];
sx q[2];
rz(0.26813771) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1902705) q[1];
sx q[1];
rz(-0.45004216) q[1];
sx q[1];
rz(-1.3954891) q[1];
x q[2];
rz(0.13161259) q[3];
sx q[3];
rz(-2.0672879) q[3];
sx q[3];
rz(2.541996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1136855) q[2];
sx q[2];
rz(-0.62782851) q[2];
sx q[2];
rz(0.17024635) q[2];
rz(-1.2804735) q[3];
sx q[3];
rz(-1.6672641) q[3];
sx q[3];
rz(1.9835651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4729507) q[0];
sx q[0];
rz(-2.1773715) q[0];
sx q[0];
rz(-2.6131795) q[0];
rz(2.2083185) q[1];
sx q[1];
rz(-0.92527881) q[1];
sx q[1];
rz(0.55567137) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2599825) q[0];
sx q[0];
rz(-2.1851808) q[0];
sx q[0];
rz(-0.23991628) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6711177) q[2];
sx q[2];
rz(-2.360376) q[2];
sx q[2];
rz(-0.90512102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56951504) q[1];
sx q[1];
rz(-1.5387427) q[1];
sx q[1];
rz(-0.7116913) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8818086) q[3];
sx q[3];
rz(-2.2385983) q[3];
sx q[3];
rz(-0.0086812191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3282503) q[2];
sx q[2];
rz(-2.3395061) q[2];
sx q[2];
rz(0.30725202) q[2];
rz(-0.79636374) q[3];
sx q[3];
rz(-2.4107404) q[3];
sx q[3];
rz(-0.097631924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.55427134) q[0];
sx q[0];
rz(-1.048943) q[0];
sx q[0];
rz(-1.3294504) q[0];
rz(0.21996552) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(-1.3185917) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585478) q[0];
sx q[0];
rz(-1.3038583) q[0];
sx q[0];
rz(0.22375317) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9302915) q[2];
sx q[2];
rz(-2.6323279) q[2];
sx q[2];
rz(2.3788521) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3941806) q[1];
sx q[1];
rz(-2.5492382) q[1];
sx q[1];
rz(2.2117531) q[1];
x q[2];
rz(-0.52826442) q[3];
sx q[3];
rz(-2.3481728) q[3];
sx q[3];
rz(-1.2906645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4465176) q[2];
sx q[2];
rz(-1.7138285) q[2];
sx q[2];
rz(-1.8565149) q[2];
rz(0.88589969) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(1.5281965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1000243) q[0];
sx q[0];
rz(-2.4055241) q[0];
sx q[0];
rz(0.31684434) q[0];
rz(2.7007804) q[1];
sx q[1];
rz(-0.79439729) q[1];
sx q[1];
rz(-2.7712834) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.646831) q[0];
sx q[0];
rz(-1.5398105) q[0];
sx q[0];
rz(1.1895183) q[0];
x q[1];
rz(-0.84443386) q[2];
sx q[2];
rz(-1.053868) q[2];
sx q[2];
rz(-0.58973613) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80659513) q[1];
sx q[1];
rz(-1.8704526) q[1];
sx q[1];
rz(0.2789168) q[1];
rz(-1.466981) q[3];
sx q[3];
rz(-1.5201638) q[3];
sx q[3];
rz(1.4736444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4050498) q[2];
sx q[2];
rz(-2.6506347) q[2];
sx q[2];
rz(1.6323818) q[2];
rz(0.80471188) q[3];
sx q[3];
rz(-1.6377621) q[3];
sx q[3];
rz(0.10778431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0263154) q[0];
sx q[0];
rz(-1.4125217) q[0];
sx q[0];
rz(-1.9052292) q[0];
rz(-2.9112877) q[1];
sx q[1];
rz(-0.31619148) q[1];
sx q[1];
rz(-2.2264623) q[1];
rz(2.438782) q[2];
sx q[2];
rz(-2.3431449) q[2];
sx q[2];
rz(0.14771067) q[2];
rz(0.95190081) q[3];
sx q[3];
rz(-1.3018082) q[3];
sx q[3];
rz(-0.34886532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
