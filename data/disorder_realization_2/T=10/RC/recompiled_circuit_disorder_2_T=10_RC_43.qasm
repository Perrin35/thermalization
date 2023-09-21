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
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(-0.50049385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252285) q[0];
sx q[0];
rz(-0.59983569) q[0];
sx q[0];
rz(-3.1347549) q[0];
x q[1];
rz(3.0481553) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(1.1364394) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11195586) q[1];
sx q[1];
rz(-1.2458548) q[1];
sx q[1];
rz(-2.3287661) q[1];
rz(1.8960564) q[3];
sx q[3];
rz(-2.6824625) q[3];
sx q[3];
rz(1.9838651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3502675) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(-1.2228489) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7704849) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(1.0789385) q[0];
rz(-1.3868015) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-0.66545495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42940419) q[0];
sx q[0];
rz(-2.5948338) q[0];
sx q[0];
rz(2.5736546) q[0];
rz(-pi) q[1];
rz(2.8198492) q[2];
sx q[2];
rz(-1.4085359) q[2];
sx q[2];
rz(0.74058796) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39311073) q[1];
sx q[1];
rz(-1.0903653) q[1];
sx q[1];
rz(-1.837681) q[1];
rz(-pi) q[2];
x q[2];
rz(0.034025107) q[3];
sx q[3];
rz(-0.861654) q[3];
sx q[3];
rz(2.2505086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5370496) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(-1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.4000777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62832075) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(0.4011641) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7064352) q[2];
sx q[2];
rz(-1.7559768) q[2];
sx q[2];
rz(2.8806825) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8850419) q[1];
sx q[1];
rz(-2.1532144) q[1];
sx q[1];
rz(-0.38862733) q[1];
x q[2];
rz(1.8758043) q[3];
sx q[3];
rz(-1.7294356) q[3];
sx q[3];
rz(-0.97914417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(-1.1714606) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(-1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(-1.6943278) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(-0.34805527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9075273) q[0];
sx q[0];
rz(-0.71279991) q[0];
sx q[0];
rz(-2.6839921) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15050998) q[2];
sx q[2];
rz(-1.9410656) q[2];
sx q[2];
rz(-0.18266695) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5665633) q[1];
sx q[1];
rz(-1.6600779) q[1];
sx q[1];
rz(0.95972285) q[1];
rz(-pi) q[2];
rz(-2.9799558) q[3];
sx q[3];
rz(-1.703754) q[3];
sx q[3];
rz(-1.1383575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7248914) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(-0.42281881) q[2];
rz(2.4041798) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(-0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2622862) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(2.2166705) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(1.8431429) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86622483) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(-2.2989681) q[0];
rz(-pi) q[1];
rz(2.5351296) q[2];
sx q[2];
rz(-1.8173373) q[2];
sx q[2];
rz(0.75070565) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7057719) q[1];
sx q[1];
rz(-2.768369) q[1];
sx q[1];
rz(0.2758287) q[1];
x q[2];
rz(2.2523746) q[3];
sx q[3];
rz(-0.71435706) q[3];
sx q[3];
rz(2.0444972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9412781) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(2.6679664) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(2.1437058) q[0];
rz(0.87431327) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(-0.46674892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40690639) q[0];
sx q[0];
rz(-0.91306251) q[0];
sx q[0];
rz(1.2975733) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42748638) q[2];
sx q[2];
rz(-1.7726521) q[2];
sx q[2];
rz(2.0919378) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0007243) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(-0.15111698) q[1];
rz(-pi) q[2];
x q[2];
rz(2.09957) q[3];
sx q[3];
rz(-0.18897945) q[3];
sx q[3];
rz(-2.6721862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(-1.5647282) q[2];
rz(-2.5148897) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(-0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108114) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683068) q[0];
sx q[0];
rz(-1.562403) q[0];
sx q[0];
rz(-0.087455672) q[0];
rz(2.1599342) q[2];
sx q[2];
rz(-2.8033211) q[2];
sx q[2];
rz(0.32064082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.067804) q[1];
sx q[1];
rz(-1.5556591) q[1];
sx q[1];
rz(3.036036) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1057304) q[3];
sx q[3];
rz(-1.8814058) q[3];
sx q[3];
rz(-0.066699337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5856813) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(-1.7377724) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(-1.9246624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4306915) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(2.8523493) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(1.8274868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1159191) q[0];
sx q[0];
rz(-2.325255) q[0];
sx q[0];
rz(-1.7220108) q[0];
rz(-pi) q[1];
rz(2.7881175) q[2];
sx q[2];
rz(-2.3620053) q[2];
sx q[2];
rz(1.0970864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60963793) q[1];
sx q[1];
rz(-1.4502118) q[1];
sx q[1];
rz(-0.054467199) q[1];
rz(-pi) q[2];
rz(-2.833509) q[3];
sx q[3];
rz(-2.3039673) q[3];
sx q[3];
rz(1.3641016) q[3];
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
rz(-2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(1.0296286) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(2.5949809) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8209936) q[0];
sx q[0];
rz(-2.5490767) q[0];
sx q[0];
rz(2.2792363) q[0];
rz(1.1999646) q[2];
sx q[2];
rz(-1.9413345) q[2];
sx q[2];
rz(-1.3818936) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0134125) q[1];
sx q[1];
rz(-1.1798522) q[1];
sx q[1];
rz(0.18545111) q[1];
rz(-pi) q[2];
rz(-2.2253753) q[3];
sx q[3];
rz(-1.1715874) q[3];
sx q[3];
rz(0.44772128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1116011) q[2];
sx q[2];
rz(-0.84838715) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(-2.4243673) q[0];
rz(-1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-2.4338914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0538941) q[0];
sx q[0];
rz(-2.2367034) q[0];
sx q[0];
rz(0.44184394) q[0];
rz(-pi) q[1];
rz(-1.3466481) q[2];
sx q[2];
rz(-1.4321616) q[2];
sx q[2];
rz(2.5773406) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.087079436) q[1];
sx q[1];
rz(-2.2594249) q[1];
sx q[1];
rz(-2.7282532) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0503065) q[3];
sx q[3];
rz(-2.4647053) q[3];
sx q[3];
rz(2.6022079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6282965) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(-2.2383402) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(-2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.306504) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(2.3256336) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(1.1076526) q[2];
sx q[2];
rz(-1.9911498) q[2];
sx q[2];
rz(-0.1291612) q[2];
rz(0.054374183) q[3];
sx q[3];
rz(-1.6118703) q[3];
sx q[3];
rz(0.76606228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];