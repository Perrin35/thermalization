OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32938862) q[0];
sx q[0];
rz(3.7731054) q[0];
sx q[0];
rz(9.4256529) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(-1.0008608) q[1];
sx q[1];
rz(0.34520087) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986818) q[0];
sx q[0];
rz(-0.094176725) q[0];
sx q[0];
rz(1.3046632) q[0];
rz(-2.7852374) q[2];
sx q[2];
rz(-0.40268597) q[2];
sx q[2];
rz(0.25549437) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6908603) q[1];
sx q[1];
rz(-1.1530515) q[1];
sx q[1];
rz(-0.91207204) q[1];
rz(1.6451938) q[3];
sx q[3];
rz(-1.5870023) q[3];
sx q[3];
rz(2.5287573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(-2.933617) q[2];
rz(-0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(-0.99288565) q[0];
rz(-1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(0.77450007) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3925288) q[0];
sx q[0];
rz(-1.3925902) q[0];
sx q[0];
rz(-3.130413) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3341137) q[2];
sx q[2];
rz(-1.1801071) q[2];
sx q[2];
rz(1.1451858) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61770536) q[1];
sx q[1];
rz(-1.3464217) q[1];
sx q[1];
rz(2.7588506) q[1];
rz(-pi) q[2];
rz(-1.1011613) q[3];
sx q[3];
rz(-1.7723933) q[3];
sx q[3];
rz(-2.2045709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.896686) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(-1.0150821) q[2];
rz(-2.8924938) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772813) q[0];
sx q[0];
rz(-2.4503777) q[0];
sx q[0];
rz(-0.15750289) q[0];
rz(1.0109673) q[1];
sx q[1];
rz(-2.8087661) q[1];
sx q[1];
rz(-3.0027622) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1942753) q[0];
sx q[0];
rz(-1.2091769) q[0];
sx q[0];
rz(0.34644923) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9857668) q[2];
sx q[2];
rz(-2.2425644) q[2];
sx q[2];
rz(2.6443554) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2646059) q[1];
sx q[1];
rz(-0.35491665) q[1];
sx q[1];
rz(2.1668129) q[1];
x q[2];
rz(-2.8582358) q[3];
sx q[3];
rz(-1.4025926) q[3];
sx q[3];
rz(0.12780549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6802754) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(-0.3581363) q[2];
rz(-2.4628468) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(-2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(0.3048234) q[0];
rz(-2.7335956) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-2.2136484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6154895) q[0];
sx q[0];
rz(-0.23357059) q[0];
sx q[0];
rz(-0.98532565) q[0];
rz(-pi) q[1];
rz(0.94522743) q[2];
sx q[2];
rz(-1.06377) q[2];
sx q[2];
rz(2.8692109) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18670652) q[1];
sx q[1];
rz(-2.2518603) q[1];
sx q[1];
rz(1.3757371) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0545066) q[3];
sx q[3];
rz(-2.1355503) q[3];
sx q[3];
rz(-1.8209344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9503595) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(-2.9534269) q[2];
rz(-1.865271) q[3];
sx q[3];
rz(-1.3151582) q[3];
sx q[3];
rz(1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(3.1222043) q[0];
rz(-2.0806606) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(-1.2394989) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1000308) q[0];
sx q[0];
rz(-1.8503555) q[0];
sx q[0];
rz(2.1651405) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2197228) q[2];
sx q[2];
rz(-0.94266329) q[2];
sx q[2];
rz(-0.79337304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8722788) q[1];
sx q[1];
rz(-1.2904823) q[1];
sx q[1];
rz(-2.6107236) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40941671) q[3];
sx q[3];
rz(-1.1077266) q[3];
sx q[3];
rz(0.45209979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1218607) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(1.3266374) q[2];
rz(2.895368) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.5354079) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(2.2703222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7428674) q[0];
sx q[0];
rz(-2.350432) q[0];
sx q[0];
rz(-1.4732811) q[0];
rz(-pi) q[1];
rz(0.75052336) q[2];
sx q[2];
rz(-2.0156246) q[2];
sx q[2];
rz(-1.9409279) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.516219) q[1];
sx q[1];
rz(-1.9587079) q[1];
sx q[1];
rz(0.078887786) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5302251) q[3];
sx q[3];
rz(-1.9690367) q[3];
sx q[3];
rz(-2.7643725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8295916) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(2.269022) q[2];
rz(-1.5729337) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0503814) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(-0.76164371) q[0];
rz(0.42539445) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(-2.2147307) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8445209) q[0];
sx q[0];
rz(-1.8064587) q[0];
sx q[0];
rz(1.8051487) q[0];
x q[1];
rz(-0.92717391) q[2];
sx q[2];
rz(-2.1774051) q[2];
sx q[2];
rz(1.7546897) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35173479) q[1];
sx q[1];
rz(-0.42966336) q[1];
sx q[1];
rz(-1.0142782) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1939932) q[3];
sx q[3];
rz(-1.5316846) q[3];
sx q[3];
rz(0.87621237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1130134) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(2.2704303) q[2];
rz(-2.8431559) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(-1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153862) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(2.723208) q[0];
rz(0.2977953) q[1];
sx q[1];
rz(-1.1957542) q[1];
sx q[1];
rz(-3.0072838) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1537841) q[0];
sx q[0];
rz(-2.1001108) q[0];
sx q[0];
rz(2.3321926) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9555126) q[2];
sx q[2];
rz(-1.969841) q[2];
sx q[2];
rz(-3.0732791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93523504) q[1];
sx q[1];
rz(-0.22769732) q[1];
sx q[1];
rz(-2.5050194) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0401434) q[3];
sx q[3];
rz(-1.3162287) q[3];
sx q[3];
rz(1.001844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(2.4995787) q[2];
rz(-1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(-1.0412019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5409656) q[0];
sx q[0];
rz(-3.0525115) q[0];
sx q[0];
rz(0.11216057) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(-2.1122011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6065858) q[0];
sx q[0];
rz(-2.2748053) q[0];
sx q[0];
rz(-0.042198472) q[0];
rz(-pi) q[1];
rz(-0.30390988) q[2];
sx q[2];
rz(-1.5329156) q[2];
sx q[2];
rz(1.5335976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9255213) q[1];
sx q[1];
rz(-0.23434429) q[1];
sx q[1];
rz(1.3804803) q[1];
rz(-2.6736027) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7182497) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(-0.038012803) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(0.14122252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24899471) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(2.7227962) q[0];
rz(-1.8839802) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(-2.9990101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6194942) q[0];
sx q[0];
rz(-1.5428679) q[0];
sx q[0];
rz(0.15837196) q[0];
rz(-2.768553) q[2];
sx q[2];
rz(-2.4786886) q[2];
sx q[2];
rz(0.66696862) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31227641) q[1];
sx q[1];
rz(-1.8354776) q[1];
sx q[1];
rz(2.1618202) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1964588) q[3];
sx q[3];
rz(-0.54900733) q[3];
sx q[3];
rz(1.6459203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7964145) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(1.4043407) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(-1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1254697) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(-3.012433) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(0.47507349) q[2];
sx q[2];
rz(-1.2550086) q[2];
sx q[2];
rz(-2.8934657) q[2];
rz(3.075243) q[3];
sx q[3];
rz(-1.8649615) q[3];
sx q[3];
rz(-2.1129114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
