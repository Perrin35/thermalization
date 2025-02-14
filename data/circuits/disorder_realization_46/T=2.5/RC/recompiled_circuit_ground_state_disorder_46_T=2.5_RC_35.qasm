OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1281066) q[0];
sx q[0];
rz(3.9689316) q[0];
sx q[0];
rz(7.6710424) q[0];
rz(0.85343051) q[1];
sx q[1];
rz(5.8813385) q[1];
sx q[1];
rz(11.452236) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6809876) q[0];
sx q[0];
rz(-1.8226242) q[0];
sx q[0];
rz(-2.2973209) q[0];
rz(0.70504712) q[2];
sx q[2];
rz(-1.7181261) q[2];
sx q[2];
rz(1.4718057) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6932363) q[1];
sx q[1];
rz(-0.020388842) q[1];
sx q[1];
rz(1.8113813) q[1];
rz(-pi) q[2];
rz(-1.2306289) q[3];
sx q[3];
rz(-2.0475884) q[3];
sx q[3];
rz(-1.7309675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9806597) q[2];
sx q[2];
rz(-2.7200343) q[2];
sx q[2];
rz(1.4483615) q[2];
rz(-0.24250044) q[3];
sx q[3];
rz(-2.226053) q[3];
sx q[3];
rz(-1.8814794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5777957) q[0];
sx q[0];
rz(-1.8255434) q[0];
sx q[0];
rz(-2.1412361) q[0];
rz(-0.73161221) q[1];
sx q[1];
rz(-1.7121406) q[1];
sx q[1];
rz(-1.9614722) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6574616) q[0];
sx q[0];
rz(-0.13741446) q[0];
sx q[0];
rz(2.4262587) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.033887788) q[2];
sx q[2];
rz(-1.8485957) q[2];
sx q[2];
rz(2.6963765) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47467642) q[1];
sx q[1];
rz(-1.70494) q[1];
sx q[1];
rz(3.139702) q[1];
rz(-pi) q[2];
rz(-0.7970771) q[3];
sx q[3];
rz(-0.38603544) q[3];
sx q[3];
rz(-0.71938709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0483027) q[2];
sx q[2];
rz(-2.2403109) q[2];
sx q[2];
rz(2.062659) q[2];
rz(-0.087873936) q[3];
sx q[3];
rz(-1.352997) q[3];
sx q[3];
rz(0.24438508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3946149) q[0];
sx q[0];
rz(-1.9786388) q[0];
sx q[0];
rz(1.8633457) q[0];
rz(-2.6784015) q[1];
sx q[1];
rz(-1.3563503) q[1];
sx q[1];
rz(2.3203704) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074479178) q[0];
sx q[0];
rz(-0.68148208) q[0];
sx q[0];
rz(1.8675748) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5089919) q[2];
sx q[2];
rz(-1.2001749) q[2];
sx q[2];
rz(-2.6573617) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48135172) q[1];
sx q[1];
rz(-1.4235745) q[1];
sx q[1];
rz(-0.46105701) q[1];
rz(0.93948313) q[3];
sx q[3];
rz(-2.2243613) q[3];
sx q[3];
rz(-2.3701649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2440196) q[2];
sx q[2];
rz(-1.579318) q[2];
sx q[2];
rz(-1.0025586) q[2];
rz(1.9133441) q[3];
sx q[3];
rz(-1.4017665) q[3];
sx q[3];
rz(1.4209411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9035852) q[0];
sx q[0];
rz(-1.0500195) q[0];
sx q[0];
rz(2.8557657) q[0];
rz(2.8803275) q[1];
sx q[1];
rz(-1.8087872) q[1];
sx q[1];
rz(-0.030698311) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26063932) q[0];
sx q[0];
rz(-2.8679823) q[0];
sx q[0];
rz(0.27283313) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2476686) q[2];
sx q[2];
rz(-1.6325762) q[2];
sx q[2];
rz(3.1368352) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49467898) q[1];
sx q[1];
rz(-2.621197) q[1];
sx q[1];
rz(1.0444276) q[1];
x q[2];
rz(0.081153374) q[3];
sx q[3];
rz(-1.8834682) q[3];
sx q[3];
rz(-2.816538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8135445) q[2];
sx q[2];
rz(-1.3406465) q[2];
sx q[2];
rz(-3.1385341) q[2];
rz(-2.3004153) q[3];
sx q[3];
rz(-2.5523461) q[3];
sx q[3];
rz(2.7868748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67741126) q[0];
sx q[0];
rz(-0.84742904) q[0];
sx q[0];
rz(-1.4671951) q[0];
rz(1.5122308) q[1];
sx q[1];
rz(-2.1758175) q[1];
sx q[1];
rz(-2.1893952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1837421) q[0];
sx q[0];
rz(-1.6410367) q[0];
sx q[0];
rz(-1.5635256) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.443794) q[2];
sx q[2];
rz(-1.2659756) q[2];
sx q[2];
rz(1.4932201) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.044301) q[1];
sx q[1];
rz(-0.11137577) q[1];
sx q[1];
rz(-1.9457413) q[1];
x q[2];
rz(-0.85785474) q[3];
sx q[3];
rz(-0.80227533) q[3];
sx q[3];
rz(2.1206534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3516922) q[2];
sx q[2];
rz(-1.60195) q[2];
sx q[2];
rz(0.44463012) q[2];
rz(0.45186684) q[3];
sx q[3];
rz(-0.63932747) q[3];
sx q[3];
rz(0.062084559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43279466) q[0];
sx q[0];
rz(-2.4169156) q[0];
sx q[0];
rz(2.2213347) q[0];
rz(2.1814749) q[1];
sx q[1];
rz(-1.3939539) q[1];
sx q[1];
rz(-2.6422909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1007784) q[0];
sx q[0];
rz(-1.7216604) q[0];
sx q[0];
rz(2.8013381) q[0];
rz(-pi) q[1];
rz(-0.22412207) q[2];
sx q[2];
rz(-0.67356985) q[2];
sx q[2];
rz(-2.3526255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0734132) q[1];
sx q[1];
rz(-1.8856032) q[1];
sx q[1];
rz(0.59609658) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4003827) q[3];
sx q[3];
rz(-1.3102207) q[3];
sx q[3];
rz(-1.8678566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7509191) q[2];
sx q[2];
rz(-1.5385188) q[2];
sx q[2];
rz(-2.1092559) q[2];
rz(-0.8647626) q[3];
sx q[3];
rz(-1.4356177) q[3];
sx q[3];
rz(2.5800956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.432935) q[0];
sx q[0];
rz(-1.189804) q[0];
sx q[0];
rz(-1.7565961) q[0];
rz(0.9224836) q[1];
sx q[1];
rz(-1.205516) q[1];
sx q[1];
rz(-0.14499697) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2534197) q[0];
sx q[0];
rz(-1.3921261) q[0];
sx q[0];
rz(-0.16797845) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85741557) q[2];
sx q[2];
rz(-1.5930297) q[2];
sx q[2];
rz(2.2043101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8783826) q[1];
sx q[1];
rz(-0.07941281) q[1];
sx q[1];
rz(1.9694049) q[1];
x q[2];
rz(2.5707158) q[3];
sx q[3];
rz(-0.8626738) q[3];
sx q[3];
rz(0.25580803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9292235) q[2];
sx q[2];
rz(-0.23304686) q[2];
sx q[2];
rz(-1.93335) q[2];
rz(0.82331795) q[3];
sx q[3];
rz(-1.9602785) q[3];
sx q[3];
rz(1.3594782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061148297) q[0];
sx q[0];
rz(-0.7875945) q[0];
sx q[0];
rz(-2.0620692) q[0];
rz(-2.1615255) q[1];
sx q[1];
rz(-2.1213687) q[1];
sx q[1];
rz(-0.10990873) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82775527) q[0];
sx q[0];
rz(-1.8640564) q[0];
sx q[0];
rz(0.37773962) q[0];
rz(-2.8632322) q[2];
sx q[2];
rz(-1.5621857) q[2];
sx q[2];
rz(-1.4573163) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4706681) q[1];
sx q[1];
rz(-1.1475855) q[1];
sx q[1];
rz(-1.2434587) q[1];
x q[2];
rz(1.0411383) q[3];
sx q[3];
rz(-1.7653437) q[3];
sx q[3];
rz(1.9499078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2073233) q[2];
sx q[2];
rz(-2.0685652) q[2];
sx q[2];
rz(0.70460021) q[2];
rz(2.8568824) q[3];
sx q[3];
rz(-2.4094818) q[3];
sx q[3];
rz(0.43340161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1339742) q[0];
sx q[0];
rz(-1.824279) q[0];
sx q[0];
rz(0.91868573) q[0];
rz(0.48565117) q[1];
sx q[1];
rz(-2.698027) q[1];
sx q[1];
rz(1.1007016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9022517) q[0];
sx q[0];
rz(-1.4318236) q[0];
sx q[0];
rz(-0.13987565) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92289779) q[2];
sx q[2];
rz(-0.37247286) q[2];
sx q[2];
rz(-0.076500741) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0079643) q[1];
sx q[1];
rz(-1.4660264) q[1];
sx q[1];
rz(1.1682577) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0555567) q[3];
sx q[3];
rz(-2.5042296) q[3];
sx q[3];
rz(-1.6364545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47621581) q[2];
sx q[2];
rz(-2.2874338) q[2];
sx q[2];
rz(2.3821135) q[2];
rz(-2.1935943) q[3];
sx q[3];
rz(-1.4576603) q[3];
sx q[3];
rz(1.52004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0976902) q[0];
sx q[0];
rz(-2.6444785) q[0];
sx q[0];
rz(2.8105766) q[0];
rz(2.379592) q[1];
sx q[1];
rz(-0.29251978) q[1];
sx q[1];
rz(-0.62823546) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2822489) q[0];
sx q[0];
rz(-2.4971136) q[0];
sx q[0];
rz(-2.6656239) q[0];
x q[1];
rz(0.64543311) q[2];
sx q[2];
rz(-2.7629768) q[2];
sx q[2];
rz(-2.9032674) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.459383) q[1];
sx q[1];
rz(-1.9193238) q[1];
sx q[1];
rz(0.63732707) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9167494) q[3];
sx q[3];
rz(-1.1802434) q[3];
sx q[3];
rz(-1.8486763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1198173) q[2];
sx q[2];
rz(-0.47406083) q[2];
sx q[2];
rz(-0.21381703) q[2];
rz(2.1553433) q[3];
sx q[3];
rz(-1.8503559) q[3];
sx q[3];
rz(-0.10828644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66961359) q[0];
sx q[0];
rz(-1.6983953) q[0];
sx q[0];
rz(-1.4629913) q[0];
rz(1.1769453) q[1];
sx q[1];
rz(-1.6564449) q[1];
sx q[1];
rz(0.0080531837) q[1];
rz(-1.0604924) q[2];
sx q[2];
rz(-2.7955187) q[2];
sx q[2];
rz(-2.1375755) q[2];
rz(1.4900653) q[3];
sx q[3];
rz(-2.2207618) q[3];
sx q[3];
rz(-0.13151463) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
