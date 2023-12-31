OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4484654) q[0];
sx q[0];
rz(-2.6187596) q[0];
sx q[0];
rz(0.62358207) q[0];
rz(3.4317598) q[1];
sx q[1];
rz(5.5640339) q[1];
sx q[1];
rz(13.066864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1815163) q[0];
sx q[0];
rz(-1.5746563) q[0];
sx q[0];
rz(0.5998248) q[0];
rz(-pi) q[1];
rz(1.7206465) q[2];
sx q[2];
rz(-1.4784001) q[2];
sx q[2];
rz(-0.42042755) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3574672) q[1];
sx q[1];
rz(-0.81175121) q[1];
sx q[1];
rz(-2.0264506) q[1];
rz(-pi) q[2];
rz(-2.0088828) q[3];
sx q[3];
rz(-1.7128908) q[3];
sx q[3];
rz(-0.70664584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7913251) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(1.9187437) q[2];
rz(-1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37110776) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(-2.0626542) q[0];
rz(-1.3868015) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(-2.4761377) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7121885) q[0];
sx q[0];
rz(-2.5948338) q[0];
sx q[0];
rz(0.56793805) q[0];
rz(-pi) q[1];
rz(1.3999248) q[2];
sx q[2];
rz(-1.8881646) q[2];
sx q[2];
rz(2.3651809) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3033777) q[1];
sx q[1];
rz(-1.3347374) q[1];
sx q[1];
rz(2.6462376) q[1];
rz(-pi) q[2];
rz(1.5311702) q[3];
sx q[3];
rz(-2.431776) q[3];
sx q[3];
rz(0.94330793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5370496) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(-1.4705307) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(-0.69141928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(0.75876045) q[0];
rz(1.8485908) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.4000777) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62832075) q[0];
sx q[0];
rz(-1.8372046) q[0];
sx q[0];
rz(2.7404286) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7064352) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(-2.8806825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6058265) q[1];
sx q[1];
rz(-1.8927791) q[1];
sx q[1];
rz(-2.1893326) q[1];
rz(-pi) q[2];
rz(-0.16617822) q[3];
sx q[3];
rz(-1.2697392) q[3];
sx q[3];
rz(-0.64134669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65163461) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(0.65762323) q[2];
rz(-1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(1.7224147) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(1.6963652) q[0];
rz(1.6943278) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(2.7935374) q[1];
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
rz(-1.939417) q[2];
sx q[2];
rz(-2.7432132) q[2];
sx q[2];
rz(0.57952651) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5665633) q[1];
sx q[1];
rz(-1.4815147) q[1];
sx q[1];
rz(-2.1818698) q[1];
rz(-pi) q[2];
rz(-0.69339852) q[3];
sx q[3];
rz(-2.9326673) q[3];
sx q[3];
rz(0.25017504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7248914) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(2.7187738) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(0.084658682) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87930644) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(2.6111531) q[0];
rz(-2.2166705) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.8431429) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4209375) q[0];
sx q[0];
rz(-1.4268095) q[0];
sx q[0];
rz(1.4074586) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2735882) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(2.1538018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0188705) q[1];
sx q[1];
rz(-1.4713305) q[1];
sx q[1];
rz(0.36032569) q[1];
rz(-pi) q[2];
rz(-2.2523746) q[3];
sx q[3];
rz(-0.71435706) q[3];
sx q[3];
rz(-2.0444972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9412781) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(2.6679664) q[2];
rz(3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(-2.2672794) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(-0.46674892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40690639) q[0];
sx q[0];
rz(-2.2285301) q[0];
sx q[0];
rz(-1.2975733) q[0];
x q[1];
rz(2.6830707) q[2];
sx q[2];
rz(-2.6715171) q[2];
sx q[2];
rz(-2.2058861) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0135865) q[1];
sx q[1];
rz(-2.9832573) q[1];
sx q[1];
rz(0.30551417) q[1];
x q[2];
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
rz(-0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5647282) q[2];
rz(2.5148897) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(-0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-0.41982857) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89824642) q[0];
sx q[0];
rz(-1.4833437) q[0];
sx q[0];
rz(1.5623708) q[0];
x q[1];
rz(1.8553472) q[2];
sx q[2];
rz(-1.7562521) q[2];
sx q[2];
rz(2.4539349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.78649) q[1];
sx q[1];
rz(-0.10663248) q[1];
sx q[1];
rz(-2.998888) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7967954) q[3];
sx q[3];
rz(-1.1296009) q[3];
sx q[3];
rz(-1.3519209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.55591136) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(-0.79706556) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(-1.9246624) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109011) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(0.28924334) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.3141059) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.700456) q[0];
sx q[0];
rz(-1.4608129) q[0];
sx q[0];
rz(0.76018795) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7881175) q[2];
sx q[2];
rz(-2.3620053) q[2];
sx q[2];
rz(-1.0970864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18409477) q[1];
sx q[1];
rz(-3.0093319) q[1];
sx q[1];
rz(1.1485419) q[1];
rz(2.833509) q[3];
sx q[3];
rz(-2.3039673) q[3];
sx q[3];
rz(-1.3641016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4079995) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(1.4286263) q[2];
rz(-0.96380487) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(-2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(-3.030581) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(0.54661173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3205991) q[0];
sx q[0];
rz(-0.59251596) q[0];
sx q[0];
rz(0.86235637) q[0];
rz(0.75066363) q[2];
sx q[2];
rz(-0.51807907) q[2];
sx q[2];
rz(-0.93875611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6555772) q[1];
sx q[1];
rz(-1.7421107) q[1];
sx q[1];
rz(1.9678712) q[1];
rz(-pi) q[2];
rz(2.6528477) q[3];
sx q[3];
rz(-2.1663323) q[3];
sx q[3];
rz(-0.83317703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(1.4477504) q[2];
rz(-1.1374121) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(-2.4243673) q[0];
rz(-1.9316797) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(0.70770121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74117888) q[0];
sx q[0];
rz(-2.3614863) q[0];
sx q[0];
rz(1.0723423) q[0];
rz(-0.14214469) q[2];
sx q[2];
rz(-1.3488349) q[2];
sx q[2];
rz(-2.1665426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51703875) q[1];
sx q[1];
rz(-2.3561764) q[1];
sx q[1];
rz(-2.024827) q[1];
rz(-pi) q[2];
rz(2.7865949) q[3];
sx q[3];
rz(-0.98155752) q[3];
sx q[3];
rz(2.0139351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5132961) q[2];
sx q[2];
rz(-2.4847023) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.306504) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(-2.03394) q[2];
sx q[2];
rz(-1.9911498) q[2];
sx q[2];
rz(-0.1291612) q[2];
rz(3.0872185) q[3];
sx q[3];
rz(-1.5297223) q[3];
sx q[3];
rz(-2.3755304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
