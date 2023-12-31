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
rz(-2.5180106) q[0];
rz(-2.8514255) q[1];
sx q[1];
rz(-0.71915141) q[1];
sx q[1];
rz(-2.6410988) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96007632) q[0];
sx q[0];
rz(-1.5746563) q[0];
sx q[0];
rz(-2.5417679) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1262769) q[2];
sx q[2];
rz(-0.17586389) q[2];
sx q[2];
rz(1.6989087) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7841255) q[1];
sx q[1];
rz(-0.81175121) q[1];
sx q[1];
rz(-1.1151421) q[1];
rz(-pi) q[2];
rz(-1.8960564) q[3];
sx q[3];
rz(-0.45913011) q[3];
sx q[3];
rz(-1.1577275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3502675) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(1.2228489) q[2];
rz(-1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(0.97035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.3868015) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(0.66545495) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9294445) q[0];
sx q[0];
rz(-1.11709) q[0];
sx q[0];
rz(-1.2544022) q[0];
x q[1];
rz(-2.6639054) q[2];
sx q[2];
rz(-0.3590695) q[2];
sx q[2];
rz(1.2815086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0004955) q[1];
sx q[1];
rz(-0.54447237) q[1];
sx q[1];
rz(2.673124) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2802248) q[3];
sx q[3];
rz(-1.5966166) q[3];
sx q[3];
rz(2.4840419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60454303) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(-3.0664505) q[2];
rz(1.6710619) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5866518) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.4000777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4979567) q[0];
sx q[0];
rz(-0.47753497) q[0];
sx q[0];
rz(-0.60995539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4351575) q[2];
sx q[2];
rz(-1.7559768) q[2];
sx q[2];
rz(2.8806825) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53576614) q[1];
sx q[1];
rz(-1.8927791) q[1];
sx q[1];
rz(-2.1893326) q[1];
rz(-1.8758043) q[3];
sx q[3];
rz(-1.7294356) q[3];
sx q[3];
rz(-2.1624485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(1.1714606) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(-1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(1.6943278) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(2.7935374) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570243) q[0];
sx q[0];
rz(-2.1977402) q[0];
sx q[0];
rz(-1.2059962) q[0];
rz(0.15050998) q[2];
sx q[2];
rz(-1.2005271) q[2];
sx q[2];
rz(-0.18266695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.066612331) q[1];
sx q[1];
rz(-0.96251026) q[1];
sx q[1];
rz(-3.0327256) q[1];
rz(-pi) q[2];
rz(0.16163687) q[3];
sx q[3];
rz(-1.703754) q[3];
sx q[3];
rz(-1.1383575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41670123) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(2.7187738) q[2];
rz(-0.73741284) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(-0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(2.6111531) q[0];
rz(0.92492217) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(-1.8431429) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4209375) q[0];
sx q[0];
rz(-1.7147831) q[0];
sx q[0];
rz(-1.7341341) q[0];
rz(-pi) q[1];
rz(-1.2735882) q[2];
sx q[2];
rz(-2.1564335) q[2];
sx q[2];
rz(2.1538018) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0188705) q[1];
sx q[1];
rz(-1.4713305) q[1];
sx q[1];
rz(0.36032569) q[1];
rz(-2.2523746) q[3];
sx q[3];
rz(-2.4272356) q[3];
sx q[3];
rz(2.0444972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(-2.6679664) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(-2.30106) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50399238) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(0.87431327) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(-0.46674892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8080374) q[0];
sx q[0];
rz(-1.7859965) q[0];
sx q[0];
rz(2.4654885) q[0];
rz(1.3495965) q[2];
sx q[2];
rz(-1.9890519) q[2];
sx q[2];
rz(-0.43005558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14086831) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(-0.15111698) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7344597) q[3];
sx q[3];
rz(-1.4758849) q[3];
sx q[3];
rz(1.519219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2816887) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(-1.5768645) q[2];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108114) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(0.41982857) q[0];
rz(-2.9201674) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5484757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89824642) q[0];
sx q[0];
rz(-1.4833437) q[0];
sx q[0];
rz(1.5792219) q[0];
rz(-pi) q[1];
rz(-0.98165841) q[2];
sx q[2];
rz(-2.8033211) q[2];
sx q[2];
rz(0.32064082) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.78649) q[1];
sx q[1];
rz(-3.0349602) q[1];
sx q[1];
rz(0.14270466) q[1];
rz(-pi) q[2];
rz(1.1057304) q[3];
sx q[3];
rz(-1.8814058) q[3];
sx q[3];
rz(3.0748933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.55591136) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(0.79706556) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4306915) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(-0.28924334) q[0];
rz(2.5166683) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.3141059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8970012) q[0];
sx q[0];
rz(-2.3750711) q[0];
sx q[0];
rz(-2.9826829) q[0];
rz(-1.9004702) q[2];
sx q[2];
rz(-2.2909819) q[2];
sx q[2];
rz(-0.61851293) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96771679) q[1];
sx q[1];
rz(-1.516725) q[1];
sx q[1];
rz(1.4500344) q[1];
rz(1.2460327) q[3];
sx q[3];
rz(-0.78403463) q[3];
sx q[3];
rz(1.3336381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(1.7129664) q[2];
rz(-0.96380487) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.2458941) q[0];
rz(-3.030581) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(2.5949809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36771691) q[0];
sx q[0];
rz(-1.1989294) q[0];
sx q[0];
rz(1.0982151) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39491744) q[2];
sx q[2];
rz(-1.2262218) q[2];
sx q[2];
rz(-0.049023703) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0134125) q[1];
sx q[1];
rz(-1.1798522) q[1];
sx q[1];
rz(-2.9561415) q[1];
rz(-pi) q[2];
rz(2.176748) q[3];
sx q[3];
rz(-2.3904739) q[3];
sx q[3];
rz(-1.5918819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(1.4477504) q[2];
rz(-1.1374121) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(-0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85957134) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(-0.71722537) q[0];
rz(-1.9316797) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(-2.4338914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087698547) q[0];
sx q[0];
rz(-2.2367034) q[0];
sx q[0];
rz(2.6997487) q[0];
rz(0.14214469) q[2];
sx q[2];
rz(-1.7927577) q[2];
sx q[2];
rz(-2.1665426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0545132) q[1];
sx q[1];
rz(-0.88216773) q[1];
sx q[1];
rz(2.7282532) q[1];
rz(-pi) q[2];
rz(-0.35499771) q[3];
sx q[3];
rz(-2.1600351) q[3];
sx q[3];
rz(-2.0139351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(-0.90325242) q[2];
rz(1.603027) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(-2.3564561) q[2];
sx q[2];
rz(-2.5265836) q[2];
sx q[2];
rz(2.1267736) q[2];
rz(-1.5296616) q[3];
sx q[3];
rz(-1.516468) q[3];
sx q[3];
rz(-0.80249912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
