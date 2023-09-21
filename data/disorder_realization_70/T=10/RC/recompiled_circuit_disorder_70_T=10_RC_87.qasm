OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6887309) q[0];
sx q[0];
rz(-2.9714669) q[0];
sx q[0];
rz(-2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2848628) q[0];
sx q[0];
rz(-0.75548178) q[0];
sx q[0];
rz(1.612624) q[0];
rz(-pi) q[1];
rz(-2.2251031) q[2];
sx q[2];
rz(-1.5208897) q[2];
sx q[2];
rz(0.33545845) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12373911) q[1];
sx q[1];
rz(-2.1542319) q[1];
sx q[1];
rz(2.6052193) q[1];
x q[2];
rz(-1.5942469) q[3];
sx q[3];
rz(-1.970885) q[3];
sx q[3];
rz(1.401702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3723345) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(-1.1713722) q[0];
rz(-0.11110335) q[2];
sx q[2];
rz(-1.9176033) q[2];
sx q[2];
rz(0.32804104) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7213388) q[1];
sx q[1];
rz(-1.4495279) q[1];
sx q[1];
rz(2.6810357) q[1];
rz(1.4731746) q[3];
sx q[3];
rz(-2.0711581) q[3];
sx q[3];
rz(2.9293159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(0.58829266) q[2];
rz(0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(0.72552848) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(0.75769889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.706447) q[0];
sx q[0];
rz(-0.62951127) q[0];
sx q[0];
rz(-2.5802617) q[0];
rz(2.4865815) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(2.8734145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1670907) q[1];
sx q[1];
rz(-2.3579512) q[1];
sx q[1];
rz(-2.2553315) q[1];
rz(-pi) q[2];
rz(0.14962872) q[3];
sx q[3];
rz(-2.4976839) q[3];
sx q[3];
rz(-0.097188918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(-2.6304723) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.3285332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851345) q[0];
sx q[0];
rz(-1.7256323) q[0];
sx q[0];
rz(-1.3550718) q[0];
rz(-pi) q[1];
rz(1.0765692) q[2];
sx q[2];
rz(-2.2708714) q[2];
sx q[2];
rz(1.9698524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7698344) q[1];
sx q[1];
rz(-1.3039939) q[1];
sx q[1];
rz(-0.00053243551) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1174326) q[3];
sx q[3];
rz(-2.8957267) q[3];
sx q[3];
rz(2.9256431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(-0.85737491) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(0.3381981) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.757526) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(1.1419883) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46576969) q[2];
sx q[2];
rz(-1.911474) q[2];
sx q[2];
rz(2.6054232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63197631) q[1];
sx q[1];
rz(-0.39034931) q[1];
sx q[1];
rz(-0.901464) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(-0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.6476691) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(-1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21022739) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(-2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-2.9528023) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81297368) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(-2.5286753) q[0];
rz(1.2994453) q[2];
sx q[2];
rz(-1.4704629) q[2];
sx q[2];
rz(0.041610418) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97555679) q[1];
sx q[1];
rz(-1.3284725) q[1];
sx q[1];
rz(1.6386599) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64340274) q[3];
sx q[3];
rz(-0.2989558) q[3];
sx q[3];
rz(1.9770196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.171689) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(-0.32399696) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(1.7623998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24810476) q[0];
sx q[0];
rz(-1.4890492) q[0];
sx q[0];
rz(-1.3260613) q[0];
x q[1];
rz(1.9786644) q[2];
sx q[2];
rz(-1.9991572) q[2];
sx q[2];
rz(-2.5025764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5768347) q[1];
sx q[1];
rz(-1.733629) q[1];
sx q[1];
rz(-0.30270438) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2077683) q[3];
sx q[3];
rz(-1.0855506) q[3];
sx q[3];
rz(-2.3683734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(-0.33106783) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(-2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11809764) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.2932628) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8338776) q[0];
sx q[0];
rz(-2.6617962) q[0];
sx q[0];
rz(-1.1307554) q[0];
rz(0.13883491) q[2];
sx q[2];
rz(-2.5359557) q[2];
sx q[2];
rz(0.56251898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5871208) q[1];
sx q[1];
rz(-1.1677824) q[1];
sx q[1];
rz(2.9906669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39956283) q[3];
sx q[3];
rz(-1.385653) q[3];
sx q[3];
rz(-3.1294587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(-0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.593489) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2244789) q[0];
sx q[0];
rz(-1.0815485) q[0];
sx q[0];
rz(-0.028562336) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31918819) q[2];
sx q[2];
rz(-0.50779283) q[2];
sx q[2];
rz(2.9920981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(-0.25823621) q[1];
rz(-pi) q[2];
rz(-2.3171114) q[3];
sx q[3];
rz(-1.5319676) q[3];
sx q[3];
rz(-2.3076434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.320497) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(2.9750032) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-2.0064158) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(0.36718711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254926) q[0];
sx q[0];
rz(-1.6187795) q[0];
sx q[0];
rz(1.539649) q[0];
rz(-2.7634002) q[2];
sx q[2];
rz(-2.4792255) q[2];
sx q[2];
rz(2.3364002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2314184) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(-0.66399666) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1807947) q[3];
sx q[3];
rz(-0.73878091) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(-1.0661351) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(0.74942855) q[2];
sx q[2];
rz(-1.5127758) q[2];
sx q[2];
rz(1.5681058) q[2];
rz(-0.84898938) q[3];
sx q[3];
rz(-1.6860387) q[3];
sx q[3];
rz(0.54308346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
