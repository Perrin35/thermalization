OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(3.9305384) q[1];
sx q[1];
rz(9.5126704) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.041931) q[0];
sx q[0];
rz(-2.4798205) q[0];
sx q[0];
rz(2.3753138) q[0];
rz(-pi) q[1];
rz(0.42638875) q[2];
sx q[2];
rz(-1.083056) q[2];
sx q[2];
rz(2.3418155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3275798) q[1];
sx q[1];
rz(-0.50427932) q[1];
sx q[1];
rz(3.1370647) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1357364) q[3];
sx q[3];
rz(-1.8294888) q[3];
sx q[3];
rz(0.70219016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(1.6128444) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(-1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.59812087) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(1.2373295) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.233261) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935826) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(2.4539024) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67458646) q[2];
sx q[2];
rz(-1.470675) q[2];
sx q[2];
rz(-1.1652201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.997628) q[1];
sx q[1];
rz(-1.3647172) q[1];
sx q[1];
rz(3.1254083) q[1];
rz(-1.3319098) q[3];
sx q[3];
rz(-0.91730648) q[3];
sx q[3];
rz(1.4373923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0040940293) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(-3.1393576) q[2];
rz(0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-2.2900443) q[0];
rz(-2.3705204) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(-0.071203701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6547346) q[0];
sx q[0];
rz(-1.808597) q[0];
sx q[0];
rz(-2.9391857) q[0];
rz(0.26239563) q[2];
sx q[2];
rz(-2.0121644) q[2];
sx q[2];
rz(2.1363832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5350266) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(1.7391298) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76020469) q[3];
sx q[3];
rz(-1.01902) q[3];
sx q[3];
rz(-2.9387568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4098542) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(-1.1666974) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.694214) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71890812) q[0];
sx q[0];
rz(-1.3410765) q[0];
sx q[0];
rz(0.16750383) q[0];
rz(1.3924613) q[2];
sx q[2];
rz(-1.5618088) q[2];
sx q[2];
rz(1.0967147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73562685) q[1];
sx q[1];
rz(-1.7298797) q[1];
sx q[1];
rz(-2.2707978) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7778266) q[3];
sx q[3];
rz(-1.7978661) q[3];
sx q[3];
rz(2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(-2.9966667) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99575627) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(2.4575535) q[0];
rz(1.1902635) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(1.2129983) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1319259) q[0];
sx q[0];
rz(-1.3103232) q[0];
sx q[0];
rz(2.7564604) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5118675) q[2];
sx q[2];
rz(-2.1246398) q[2];
sx q[2];
rz(0.50309203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4124406) q[1];
sx q[1];
rz(-1.3153331) q[1];
sx q[1];
rz(1.4490119) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8323684) q[3];
sx q[3];
rz(-1.7422758) q[3];
sx q[3];
rz(-1.24303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6325536) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(-0.63009134) q[2];
rz(-2.5693494) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(-0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.1155788) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(-0.28636006) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(0.58475959) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2858546) q[0];
sx q[0];
rz(-1.2585088) q[0];
sx q[0];
rz(0.62494846) q[0];
x q[1];
rz(-0.61442394) q[2];
sx q[2];
rz(-1.5493869) q[2];
sx q[2];
rz(1.1504088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44536351) q[1];
sx q[1];
rz(-0.74659691) q[1];
sx q[1];
rz(2.0018875) q[1];
x q[2];
rz(-1.6213958) q[3];
sx q[3];
rz(-1.4046602) q[3];
sx q[3];
rz(-1.43464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(-0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0528089) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(2.5928296) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(2.5351977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8618146) q[0];
sx q[0];
rz(-1.8068411) q[0];
sx q[0];
rz(-1.0311014) q[0];
rz(-0.90784448) q[2];
sx q[2];
rz(-1.8722476) q[2];
sx q[2];
rz(-2.9482258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8555774) q[1];
sx q[1];
rz(-2.3708323) q[1];
sx q[1];
rz(-1.7325749) q[1];
rz(-pi) q[2];
rz(1.1612732) q[3];
sx q[3];
rz(-1.8990574) q[3];
sx q[3];
rz(1.8144153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(-2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33928076) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(0.75396496) q[0];
rz(-2.1915961) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(-0.22769134) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2988417) q[0];
sx q[0];
rz(-1.9831053) q[0];
sx q[0];
rz(3.0639248) q[0];
rz(-pi) q[1];
rz(1.20485) q[2];
sx q[2];
rz(-1.9908675) q[2];
sx q[2];
rz(-2.7811188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88340064) q[1];
sx q[1];
rz(-1.0419783) q[1];
sx q[1];
rz(-1.3099111) q[1];
x q[2];
rz(-2.5826107) q[3];
sx q[3];
rz(-1.535927) q[3];
sx q[3];
rz(2.1502286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(2.3030346) q[2];
rz(-2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39695981) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(0.80378419) q[0];
rz(2.0869758) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(1.9931591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20439416) q[0];
sx q[0];
rz(-0.075852595) q[0];
sx q[0];
rz(1.2501459) q[0];
rz(-pi) q[1];
rz(-2.7788413) q[2];
sx q[2];
rz(-2.0667549) q[2];
sx q[2];
rz(2.9438058) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25363898) q[1];
sx q[1];
rz(-1.9983074) q[1];
sx q[1];
rz(-1.0484139) q[1];
rz(-pi) q[2];
rz(-0.095454772) q[3];
sx q[3];
rz(-2.7881552) q[3];
sx q[3];
rz(-0.27775882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-2.3633374) q[2];
rz(-2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2739928) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(-0.25319779) q[0];
rz(2.4018438) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(-2.443312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90726501) q[0];
sx q[0];
rz(-2.8303879) q[0];
sx q[0];
rz(-0.23647232) q[0];
x q[1];
rz(-1.8928705) q[2];
sx q[2];
rz(-0.078157166) q[2];
sx q[2];
rz(0.31101481) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5729534) q[1];
sx q[1];
rz(-1.1146953) q[1];
sx q[1];
rz(-2.9005463) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0032759) q[3];
sx q[3];
rz(-1.7756117) q[3];
sx q[3];
rz(1.8703465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0239821) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(-2.3013766) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(-1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8828076) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(3.1124658) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(1.2226979) q[2];
sx q[2];
rz(-0.99708996) q[2];
sx q[2];
rz(-0.70380824) q[2];
rz(-1.7520262) q[3];
sx q[3];
rz(-2.443234) q[3];
sx q[3];
rz(0.20886226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
