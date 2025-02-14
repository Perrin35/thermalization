OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4360566) q[0];
sx q[0];
rz(5.4242791) q[0];
sx q[0];
rz(10.407596) q[0];
rz(-1.9614027) q[1];
sx q[1];
rz(-0.72620121) q[1];
sx q[1];
rz(0.89816165) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9010278) q[0];
sx q[0];
rz(-1.8870346) q[0];
sx q[0];
rz(-1.4610806) q[0];
rz(-pi) q[1];
rz(-2.5835393) q[2];
sx q[2];
rz(-0.61667569) q[2];
sx q[2];
rz(1.0215173) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4136988) q[1];
sx q[1];
rz(-2.4174712) q[1];
sx q[1];
rz(-0.65774386) q[1];
rz(-1.6706561) q[3];
sx q[3];
rz(-1.4263527) q[3];
sx q[3];
rz(2.5564872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0183705) q[2];
sx q[2];
rz(-1.3192588) q[2];
sx q[2];
rz(2.3940864) q[2];
rz(-1.2627164) q[3];
sx q[3];
rz(-1.5371753) q[3];
sx q[3];
rz(0.023716299) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30775192) q[0];
sx q[0];
rz(-3.096014) q[0];
sx q[0];
rz(1.0652834) q[0];
rz(-2.4899958) q[1];
sx q[1];
rz(-0.35531303) q[1];
sx q[1];
rz(-1.5078872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9621443) q[0];
sx q[0];
rz(-2.1198258) q[0];
sx q[0];
rz(0.35513504) q[0];
rz(-pi) q[1];
rz(1.2149732) q[2];
sx q[2];
rz(-1.1220555) q[2];
sx q[2];
rz(-0.85576487) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7439241) q[1];
sx q[1];
rz(-0.63536352) q[1];
sx q[1];
rz(-0.85224908) q[1];
rz(-1.5298858) q[3];
sx q[3];
rz(-2.5299151) q[3];
sx q[3];
rz(3.1168695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4773341) q[2];
sx q[2];
rz(-1.6691672) q[2];
sx q[2];
rz(3.0490457) q[2];
rz(0.44068286) q[3];
sx q[3];
rz(-0.40658545) q[3];
sx q[3];
rz(1.996076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6630845) q[0];
sx q[0];
rz(-0.19190754) q[0];
sx q[0];
rz(-1.2364016) q[0];
rz(2.9036486) q[1];
sx q[1];
rz(-1.4711719) q[1];
sx q[1];
rz(2.1740289) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9151354) q[0];
sx q[0];
rz(-0.98513705) q[0];
sx q[0];
rz(2.4047635) q[0];
rz(2.6009623) q[2];
sx q[2];
rz(-1.2714567) q[2];
sx q[2];
rz(0.23330748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1246207) q[1];
sx q[1];
rz(-1.3836814) q[1];
sx q[1];
rz(-1.2618503) q[1];
rz(-1.1573575) q[3];
sx q[3];
rz(-1.6235768) q[3];
sx q[3];
rz(2.6708713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48270109) q[2];
sx q[2];
rz(-0.37223688) q[2];
sx q[2];
rz(2.2230395) q[2];
rz(0.77361584) q[3];
sx q[3];
rz(-2.2548803) q[3];
sx q[3];
rz(-0.95019597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8755662) q[0];
sx q[0];
rz(-0.73930621) q[0];
sx q[0];
rz(-2.8493122) q[0];
rz(0.87458163) q[1];
sx q[1];
rz(-0.62868172) q[1];
sx q[1];
rz(0.94327092) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1312201) q[0];
sx q[0];
rz(-1.2772075) q[0];
sx q[0];
rz(-2.2468727) q[0];
x q[1];
rz(1.1601527) q[2];
sx q[2];
rz(-2.6371397) q[2];
sx q[2];
rz(-1.8093833) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29674655) q[1];
sx q[1];
rz(-1.2592717) q[1];
sx q[1];
rz(-0.89235899) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1481879) q[3];
sx q[3];
rz(-1.3311989) q[3];
sx q[3];
rz(1.8519608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6167407) q[2];
sx q[2];
rz(-2.0276232) q[2];
sx q[2];
rz(-2.2576766) q[2];
rz(-1.135745) q[3];
sx q[3];
rz(-0.92848778) q[3];
sx q[3];
rz(2.4511724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8166872) q[0];
sx q[0];
rz(-0.08520928) q[0];
sx q[0];
rz(0.19164044) q[0];
rz(-1.8457671) q[1];
sx q[1];
rz(-1.192966) q[1];
sx q[1];
rz(0.4506909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86203456) q[0];
sx q[0];
rz(-2.4992895) q[0];
sx q[0];
rz(-1.7234283) q[0];
rz(-pi) q[1];
rz(-0.27235548) q[2];
sx q[2];
rz(-1.5030454) q[2];
sx q[2];
rz(2.1960047) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60495341) q[1];
sx q[1];
rz(-0.99483314) q[1];
sx q[1];
rz(0.085809068) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21416625) q[3];
sx q[3];
rz(-2.8422822) q[3];
sx q[3];
rz(2.4714937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6777163) q[2];
sx q[2];
rz(-1.8666942) q[2];
sx q[2];
rz(0.95787588) q[2];
rz(3.0887582) q[3];
sx q[3];
rz(-0.54532471) q[3];
sx q[3];
rz(2.7132645) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65474725) q[0];
sx q[0];
rz(-1.0423648) q[0];
sx q[0];
rz(2.6690707) q[0];
rz(-0.76332244) q[1];
sx q[1];
rz(-2.4210052) q[1];
sx q[1];
rz(-0.29019132) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5369072) q[0];
sx q[0];
rz(-1.3504656) q[0];
sx q[0];
rz(-2.4308886) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1212513) q[2];
sx q[2];
rz(-1.2977306) q[2];
sx q[2];
rz(-2.9228763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7136894) q[1];
sx q[1];
rz(-0.91304243) q[1];
sx q[1];
rz(1.7005928) q[1];
rz(-pi) q[2];
rz(-1.7700341) q[3];
sx q[3];
rz(-1.3955978) q[3];
sx q[3];
rz(-0.74149473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.36576888) q[2];
sx q[2];
rz(-0.86161986) q[2];
sx q[2];
rz(-2.5422868) q[2];
rz(-1.722909) q[3];
sx q[3];
rz(-1.8214046) q[3];
sx q[3];
rz(2.1883709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80239427) q[0];
sx q[0];
rz(-1.1856439) q[0];
sx q[0];
rz(-1.8307357) q[0];
rz(1.6400379) q[1];
sx q[1];
rz(-2.7412667) q[1];
sx q[1];
rz(-0.60360533) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5854322) q[0];
sx q[0];
rz(-1.4370262) q[0];
sx q[0];
rz(1.9942392) q[0];
x q[1];
rz(1.3519751) q[2];
sx q[2];
rz(-1.0299152) q[2];
sx q[2];
rz(2.4336124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7052066) q[1];
sx q[1];
rz(-1.6387741) q[1];
sx q[1];
rz(1.3297362) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9333618) q[3];
sx q[3];
rz(-2.3891267) q[3];
sx q[3];
rz(-0.83112874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32460585) q[2];
sx q[2];
rz(-0.88674712) q[2];
sx q[2];
rz(2.9465607) q[2];
rz(2.4845691) q[3];
sx q[3];
rz(-1.7200836) q[3];
sx q[3];
rz(3.0278964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8103771) q[0];
sx q[0];
rz(-1.4137784) q[0];
sx q[0];
rz(3.069416) q[0];
rz(2.3585034) q[1];
sx q[1];
rz(-0.63242811) q[1];
sx q[1];
rz(-0.91792387) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4179787) q[0];
sx q[0];
rz(-2.0416913) q[0];
sx q[0];
rz(1.3379452) q[0];
rz(-pi) q[1];
rz(0.36717271) q[2];
sx q[2];
rz(-2.4879666) q[2];
sx q[2];
rz(3.0116561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.16282475) q[1];
sx q[1];
rz(-0.37436327) q[1];
sx q[1];
rz(2.4434213) q[1];
rz(-1.9959171) q[3];
sx q[3];
rz(-1.8048058) q[3];
sx q[3];
rz(-2.0834415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5607295) q[2];
sx q[2];
rz(-1.9274638) q[2];
sx q[2];
rz(-1.5011935) q[2];
rz(-0.8864657) q[3];
sx q[3];
rz(-1.3366046) q[3];
sx q[3];
rz(0.10543536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5378872) q[0];
sx q[0];
rz(-0.34905809) q[0];
sx q[0];
rz(2.6339997) q[0];
rz(-0.72788584) q[1];
sx q[1];
rz(-1.6517086) q[1];
sx q[1];
rz(-1.1061888) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87306753) q[0];
sx q[0];
rz(-1.2374479) q[0];
sx q[0];
rz(0.26645904) q[0];
x q[1];
rz(2.699643) q[2];
sx q[2];
rz(-0.39296752) q[2];
sx q[2];
rz(3.0073187) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0303846) q[1];
sx q[1];
rz(-1.9002575) q[1];
sx q[1];
rz(-1.9558286) q[1];
x q[2];
rz(0.9000128) q[3];
sx q[3];
rz(-2.4379895) q[3];
sx q[3];
rz(2.7629064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20405208) q[2];
sx q[2];
rz(-1.1037408) q[2];
sx q[2];
rz(-0.44720116) q[2];
rz(1.7043381) q[3];
sx q[3];
rz(-2.3279326) q[3];
sx q[3];
rz(1.6566431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60780418) q[0];
sx q[0];
rz(-2.0543126) q[0];
sx q[0];
rz(2.6265889) q[0];
rz(1.1732514) q[1];
sx q[1];
rz(-1.8517905) q[1];
sx q[1];
rz(0.44949964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363578) q[0];
sx q[0];
rz(-0.32669386) q[0];
sx q[0];
rz(1.5004083) q[0];
rz(-pi) q[1];
rz(2.1859288) q[2];
sx q[2];
rz(-2.0086096) q[2];
sx q[2];
rz(2.4085338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0970407) q[1];
sx q[1];
rz(-2.1637193) q[1];
sx q[1];
rz(-2.135599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24522606) q[3];
sx q[3];
rz(-0.86638993) q[3];
sx q[3];
rz(-2.3622262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0476734) q[2];
sx q[2];
rz(-2.8670222) q[2];
sx q[2];
rz(0.74335113) q[2];
rz(-0.36176935) q[3];
sx q[3];
rz(-2.1105364) q[3];
sx q[3];
rz(0.37775347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34306985) q[0];
sx q[0];
rz(-1.1953851) q[0];
sx q[0];
rz(0.60085798) q[0];
rz(-2.9877904) q[1];
sx q[1];
rz(-1.8228795) q[1];
sx q[1];
rz(0.22620329) q[1];
rz(-1.9696196) q[2];
sx q[2];
rz(-2.8097967) q[2];
sx q[2];
rz(1.9423165) q[2];
rz(-1.947851) q[3];
sx q[3];
rz(-0.98529639) q[3];
sx q[3];
rz(-2.1724971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
