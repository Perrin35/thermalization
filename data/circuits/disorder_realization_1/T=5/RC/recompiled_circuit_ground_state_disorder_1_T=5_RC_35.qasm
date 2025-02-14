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
rz(0.57408339) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.047217) q[0];
sx q[0];
rz(-1.2021016) q[0];
sx q[0];
rz(-2.019702) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98784222) q[2];
sx q[2];
rz(-1.541809) q[2];
sx q[2];
rz(-1.0430973) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8694386) q[1];
sx q[1];
rz(-1.4813786) q[1];
sx q[1];
rz(-1.2050259) q[1];
rz(-pi) q[2];
rz(-1.194379) q[3];
sx q[3];
rz(-1.5113081) q[3];
sx q[3];
rz(0.74947442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2744039) q[2];
sx q[2];
rz(-1.8827266) q[2];
sx q[2];
rz(-1.3192419) q[2];
rz(0.073171767) q[3];
sx q[3];
rz(-1.3061482) q[3];
sx q[3];
rz(-0.81899548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0935593) q[0];
sx q[0];
rz(-1.9608542) q[0];
sx q[0];
rz(2.4617885) q[0];
rz(-0.80348429) q[1];
sx q[1];
rz(-1.7796703) q[1];
sx q[1];
rz(0.63978535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073994324) q[0];
sx q[0];
rz(-0.72866458) q[0];
sx q[0];
rz(2.0557899) q[0];
rz(-pi) q[1];
rz(-0.72987188) q[2];
sx q[2];
rz(-1.6449071) q[2];
sx q[2];
rz(-0.48219901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8315994) q[1];
sx q[1];
rz(-0.91504708) q[1];
sx q[1];
rz(0.80226957) q[1];
x q[2];
rz(-1.3469668) q[3];
sx q[3];
rz(-1.8719684) q[3];
sx q[3];
rz(-0.45066842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5378319) q[2];
sx q[2];
rz(-0.49850264) q[2];
sx q[2];
rz(2.4375088) q[2];
rz(-3.0294561) q[3];
sx q[3];
rz(-1.6928558) q[3];
sx q[3];
rz(2.1191547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0391487) q[0];
sx q[0];
rz(-1.1071438) q[0];
sx q[0];
rz(2.802134) q[0];
rz(2.0231694) q[1];
sx q[1];
rz(-2.2148841) q[1];
sx q[1];
rz(-0.035645398) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4887071) q[0];
sx q[0];
rz(-1.4849326) q[0];
sx q[0];
rz(2.9339132) q[0];
rz(-0.33907922) q[2];
sx q[2];
rz(-1.2847752) q[2];
sx q[2];
rz(1.3468483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.798637) q[1];
sx q[1];
rz(-1.3561212) q[1];
sx q[1];
rz(0.66185419) q[1];
x q[2];
rz(-2.527329) q[3];
sx q[3];
rz(-1.8395367) q[3];
sx q[3];
rz(-2.4128071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3493335) q[2];
sx q[2];
rz(-2.4534241) q[2];
sx q[2];
rz(-0.86000693) q[2];
rz(0.79259121) q[3];
sx q[3];
rz(-2.9383797) q[3];
sx q[3];
rz(-2.121076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0122796) q[0];
sx q[0];
rz(-1.3849994) q[0];
sx q[0];
rz(2.5087575) q[0];
rz(-1.9889779) q[1];
sx q[1];
rz(-2.2668138) q[1];
sx q[1];
rz(-1.4702183) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33575422) q[0];
sx q[0];
rz(-2.4921761) q[0];
sx q[0];
rz(-1.5494256) q[0];
rz(2.987902) q[2];
sx q[2];
rz(-2.3665303) q[2];
sx q[2];
rz(-2.3714921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8006261) q[1];
sx q[1];
rz(-1.8273841) q[1];
sx q[1];
rz(-0.8534732) q[1];
rz(-pi) q[2];
rz(-2.4012864) q[3];
sx q[3];
rz(-2.3362219) q[3];
sx q[3];
rz(0.80163085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4319438) q[2];
sx q[2];
rz(-2.180763) q[2];
sx q[2];
rz(2.7107837) q[2];
rz(-0.68239799) q[3];
sx q[3];
rz(-1.4902481) q[3];
sx q[3];
rz(-0.94990134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7948941) q[0];
sx q[0];
rz(-1.8759202) q[0];
sx q[0];
rz(-1.1899813) q[0];
rz(-2.8311912) q[1];
sx q[1];
rz(-1.0810532) q[1];
sx q[1];
rz(-2.6365872) q[1];
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
x q[1];
rz(-0.1311028) q[2];
sx q[2];
rz(-1.982568) q[2];
sx q[2];
rz(0.51136651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0934) q[1];
sx q[1];
rz(-1.6580771) q[1];
sx q[1];
rz(0.49836664) q[1];
rz(-3.04516) q[3];
sx q[3];
rz(-2.4135906) q[3];
sx q[3];
rz(-1.6681287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7348822) q[2];
sx q[2];
rz(-0.70047417) q[2];
sx q[2];
rz(0.90275466) q[2];
rz(1.0550176) q[3];
sx q[3];
rz(-0.86692923) q[3];
sx q[3];
rz(-2.6796851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2257253) q[0];
sx q[0];
rz(-0.71284717) q[0];
sx q[0];
rz(2.0676887) q[0];
rz(1.7621) q[1];
sx q[1];
rz(-0.72934377) q[1];
sx q[1];
rz(-3.0440547) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2052692) q[0];
sx q[0];
rz(-0.57010929) q[0];
sx q[0];
rz(-2.8015215) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1937134) q[2];
sx q[2];
rz(-1.4137739) q[2];
sx q[2];
rz(-2.825001) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3055468) q[1];
sx q[1];
rz(-1.6191481) q[1];
sx q[1];
rz(-1.6319066) q[1];
rz(-pi) q[2];
rz(2.318368) q[3];
sx q[3];
rz(-2.4308017) q[3];
sx q[3];
rz(-1.8985793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0292616) q[2];
sx q[2];
rz(-1.6882201) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77808648) q[0];
sx q[0];
rz(-0.030487617) q[0];
sx q[0];
rz(-1.7331069) q[0];
rz(2.1259437) q[1];
sx q[1];
rz(-1.8135095) q[1];
sx q[1];
rz(1.4782864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625945) q[0];
sx q[0];
rz(-1.2752879) q[0];
sx q[0];
rz(2.5062755) q[0];
x q[1];
rz(-2.4934216) q[2];
sx q[2];
rz(-1.0888466) q[2];
sx q[2];
rz(-2.3881557) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95132213) q[1];
sx q[1];
rz(-0.45004216) q[1];
sx q[1];
rz(1.3954891) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8084548) q[3];
sx q[3];
rz(-2.6293652) q[3];
sx q[3];
rz(0.32853261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1136855) q[2];
sx q[2];
rz(-0.62782851) q[2];
sx q[2];
rz(0.17024635) q[2];
rz(1.2804735) q[3];
sx q[3];
rz(-1.6672641) q[3];
sx q[3];
rz(-1.9835651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4729507) q[0];
sx q[0];
rz(-2.1773715) q[0];
sx q[0];
rz(2.6131795) q[0];
rz(-0.93327418) q[1];
sx q[1];
rz(-2.2163138) q[1];
sx q[1];
rz(2.5859213) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9708722) q[0];
sx q[0];
rz(-1.375388) q[0];
sx q[0];
rz(0.94265818) q[0];
rz(-pi) q[1];
rz(-2.3494928) q[2];
sx q[2];
rz(-1.5002155) q[2];
sx q[2];
rz(2.5472699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.96413999) q[1];
sx q[1];
rz(-0.71228668) q[1];
sx q[1];
rz(0.049055462) q[1];
rz(1.2597841) q[3];
sx q[3];
rz(-0.90299435) q[3];
sx q[3];
rz(-0.0086812191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8133424) q[2];
sx q[2];
rz(-2.3395061) q[2];
sx q[2];
rz(-2.8343406) q[2];
rz(-0.79636374) q[3];
sx q[3];
rz(-2.4107404) q[3];
sx q[3];
rz(-0.097631924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55427134) q[0];
sx q[0];
rz(-2.0926496) q[0];
sx q[0];
rz(-1.8121423) q[0];
rz(2.9216271) q[1];
sx q[1];
rz(-0.34383067) q[1];
sx q[1];
rz(1.8230009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7466965) q[0];
sx q[0];
rz(-1.7864972) q[0];
sx q[0];
rz(-1.8442276) q[0];
rz(-0.49974738) q[2];
sx q[2];
rz(-1.4683654) q[2];
sx q[2];
rz(-0.62291716) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66178759) q[1];
sx q[1];
rz(-1.1068245) q[1];
sx q[1];
rz(0.3825835) q[1];
x q[2];
rz(-1.0974466) q[3];
sx q[3];
rz(-0.90765491) q[3];
sx q[3];
rz(-1.1569661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4465176) q[2];
sx q[2];
rz(-1.7138285) q[2];
sx q[2];
rz(1.2850777) q[2];
rz(2.255693) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(1.6133962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0415683) q[0];
sx q[0];
rz(-0.73606857) q[0];
sx q[0];
rz(-0.31684434) q[0];
rz(2.7007804) q[1];
sx q[1];
rz(-2.3471954) q[1];
sx q[1];
rz(2.7712834) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.646831) q[0];
sx q[0];
rz(-1.5398105) q[0];
sx q[0];
rz(-1.9520743) q[0];
rz(0.86285186) q[2];
sx q[2];
rz(-2.2782695) q[2];
sx q[2];
rz(-1.6528636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1052484) q[1];
sx q[1];
rz(-2.7351028) q[1];
sx q[1];
rz(0.84291934) q[1];
x q[2];
rz(-2.0256151) q[3];
sx q[3];
rz(-3.0261281) q[3];
sx q[3];
rz(-2.786557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.73654282) q[2];
sx q[2];
rz(-2.6506347) q[2];
sx q[2];
rz(1.5092108) q[2];
rz(-2.3368808) q[3];
sx q[3];
rz(-1.6377621) q[3];
sx q[3];
rz(-3.0338083) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1152773) q[0];
sx q[0];
rz(-1.4125217) q[0];
sx q[0];
rz(-1.9052292) q[0];
rz(-0.23030494) q[1];
sx q[1];
rz(-2.8254012) q[1];
sx q[1];
rz(0.91513035) q[1];
rz(-0.7028107) q[2];
sx q[2];
rz(-2.3431449) q[2];
sx q[2];
rz(0.14771067) q[2];
rz(-2.8152499) q[3];
sx q[3];
rz(-0.97728609) q[3];
sx q[3];
rz(1.4090007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
