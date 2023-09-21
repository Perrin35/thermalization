OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(-0.52283302) q[0];
sx q[0];
rz(-0.62358207) q[0];
rz(3.4317598) q[1];
sx q[1];
rz(5.5640339) q[1];
sx q[1];
rz(13.066864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1815163) q[0];
sx q[0];
rz(-1.5669364) q[0];
sx q[0];
rz(0.5998248) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0481553) q[2];
sx q[2];
rz(-1.7200025) q[2];
sx q[2];
rz(-2.0051533) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9762293) q[1];
sx q[1];
rz(-0.86127087) q[1];
sx q[1];
rz(2.707259) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1327098) q[3];
sx q[3];
rz(-1.7128908) q[3];
sx q[3];
rz(-0.70664584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7913251) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(1.2228489) q[2];
rz(1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(2.1712415) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704849) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(-2.0626542) q[0];
rz(1.7547912) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-0.66545495) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21214813) q[0];
sx q[0];
rz(-1.11709) q[0];
sx q[0];
rz(1.2544022) q[0];
x q[1];
rz(-2.8198492) q[2];
sx q[2];
rz(-1.7330568) q[2];
sx q[2];
rz(0.74058796) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14109719) q[1];
sx q[1];
rz(-2.5971203) q[1];
sx q[1];
rz(-2.673124) q[1];
x q[2];
rz(-1.5311702) q[3];
sx q[3];
rz(-2.431776) q[3];
sx q[3];
rz(-0.94330793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60454303) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(3.0664505) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(-1.8485908) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.4000777) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62832075) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(-2.7404286) q[0];
rz(-pi) q[1];
rz(-2.954735) q[2];
sx q[2];
rz(-1.7041022) q[2];
sx q[2];
rz(-1.8568298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53576614) q[1];
sx q[1];
rz(-1.8927791) q[1];
sx q[1];
rz(-0.95226007) q[1];
rz(-0.16617822) q[3];
sx q[3];
rz(-1.2697392) q[3];
sx q[3];
rz(-0.64134669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(-1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(-1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(1.6963652) q[0];
rz(-1.6943278) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(-0.34805527) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30663438) q[0];
sx q[0];
rz(-1.277703) q[0];
sx q[0];
rz(-0.65960633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9449171) q[2];
sx q[2];
rz(-1.4305563) q[2];
sx q[2];
rz(-1.3333048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0749803) q[1];
sx q[1];
rz(-2.1790824) q[1];
sx q[1];
rz(3.0327256) q[1];
rz(-pi) q[2];
rz(1.4361037) q[3];
sx q[3];
rz(-1.7309942) q[3];
sx q[3];
rz(-2.687541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7248914) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(0.42281881) q[2];
rz(2.4041798) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(-0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87930644) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(2.6111531) q[0];
rz(2.2166705) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(-1.8431429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7206551) q[0];
sx q[0];
rz(-1.7147831) q[0];
sx q[0];
rz(1.4074586) q[0];
rz(-1.2735882) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(0.98779087) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1227222) q[1];
sx q[1];
rz(-1.6702622) q[1];
sx q[1];
rz(0.36032569) q[1];
x q[2];
rz(-2.2523746) q[3];
sx q[3];
rz(-2.4272356) q[3];
sx q[3];
rz(-1.0970955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(0.099362699) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(0.84053269) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50399238) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(-2.1437058) q[0];
rz(-0.87431327) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(0.46674892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8080374) q[0];
sx q[0];
rz(-1.7859965) q[0];
sx q[0];
rz(-0.67610418) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3495965) q[2];
sx q[2];
rz(-1.1525407) q[2];
sx q[2];
rz(0.43005558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0135865) q[1];
sx q[1];
rz(-2.9832573) q[1];
sx q[1];
rz(2.8360785) q[1];
rz(-pi) q[2];
rz(-2.09957) q[3];
sx q[3];
rz(-0.18897945) q[3];
sx q[3];
rz(2.6721862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5768645) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(2.4600162) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(2.9201674) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89824642) q[0];
sx q[0];
rz(-1.6582489) q[0];
sx q[0];
rz(1.5623708) q[0];
rz(-pi) q[1];
rz(-0.98165841) q[2];
sx q[2];
rz(-0.33827153) q[2];
sx q[2];
rz(-0.32064082) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6429813) q[1];
sx q[1];
rz(-1.6763408) q[1];
sx q[1];
rz(1.5860182) q[1];
rz(-pi) q[2];
rz(0.34479721) q[3];
sx q[3];
rz(-1.1296009) q[3];
sx q[3];
rz(1.7896717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55591136) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(-0.79706556) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.9246624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109011) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(-0.28924334) q[0];
rz(-2.5166683) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.8274868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1159191) q[0];
sx q[0];
rz(-2.325255) q[0];
sx q[0];
rz(-1.4195819) q[0];
x q[1];
rz(0.35347519) q[2];
sx q[2];
rz(-2.3620053) q[2];
sx q[2];
rz(-1.0970864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.96771679) q[1];
sx q[1];
rz(-1.516725) q[1];
sx q[1];
rz(-1.6915583) q[1];
rz(1.8955599) q[3];
sx q[3];
rz(-0.78403463) q[3];
sx q[3];
rz(-1.3336381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4079995) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.4286263) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52255094) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(1.2458941) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(-0.54661173) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3205991) q[0];
sx q[0];
rz(-0.59251596) q[0];
sx q[0];
rz(-2.2792363) q[0];
x q[1];
rz(-2.390929) q[2];
sx q[2];
rz(-2.6235136) q[2];
sx q[2];
rz(-2.2028365) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1281801) q[1];
sx q[1];
rz(-1.9617404) q[1];
sx q[1];
rz(-2.9561415) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96484465) q[3];
sx q[3];
rz(-2.3904739) q[3];
sx q[3];
rz(1.5497108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(-1.4477504) q[2];
rz(-1.1374121) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(2.494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820213) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(2.4243673) q[0];
rz(1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-0.70770121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1987863) q[0];
sx q[0];
rz(-1.9137303) q[0];
sx q[0];
rz(0.85533157) q[0];
x q[1];
rz(-1.7949445) q[2];
sx q[2];
rz(-1.7094311) q[2];
sx q[2];
rz(2.5773406) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0545132) q[1];
sx q[1];
rz(-2.2594249) q[1];
sx q[1];
rz(-2.7282532) q[1];
x q[2];
rz(1.0912861) q[3];
sx q[3];
rz(-2.4647053) q[3];
sx q[3];
rz(-0.53938473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6282965) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(2.2383402) q[2];
rz(1.5385657) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(-0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(-0.8159591) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(-1.1076526) q[2];
sx q[2];
rz(-1.1504428) q[2];
sx q[2];
rz(3.0124315) q[2];
rz(-0.64745263) q[3];
sx q[3];
rz(-0.068131937) q[3];
sx q[3];
rz(1.6905231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
