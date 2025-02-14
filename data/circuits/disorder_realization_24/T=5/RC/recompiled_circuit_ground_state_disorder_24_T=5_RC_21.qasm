OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4234023) q[0];
sx q[0];
rz(-0.0027522491) q[0];
sx q[0];
rz(2.1299074) q[0];
rz(0.45971316) q[1];
sx q[1];
rz(-1.8399532) q[1];
sx q[1];
rz(-0.17170061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7778439) q[0];
sx q[0];
rz(-0.96723377) q[0];
sx q[0];
rz(1.6700755) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9471875) q[2];
sx q[2];
rz(-1.462359) q[2];
sx q[2];
rz(-0.51314236) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5039035) q[1];
sx q[1];
rz(-1.1451964) q[1];
sx q[1];
rz(-0.63987672) q[1];
rz(1.0461476) q[3];
sx q[3];
rz(-3.0339315) q[3];
sx q[3];
rz(2.8058509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95919886) q[2];
sx q[2];
rz(-2.6540519) q[2];
sx q[2];
rz(-1.4200776) q[2];
rz(-1.1848263) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083953388) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(-0.49348304) q[0];
rz(-0.77589846) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(-0.50481558) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1349943) q[0];
sx q[0];
rz(-0.93824832) q[0];
sx q[0];
rz(0.68944864) q[0];
rz(-1.8680598) q[2];
sx q[2];
rz(-1.0631764) q[2];
sx q[2];
rz(1.0647237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82693726) q[1];
sx q[1];
rz(-0.75410226) q[1];
sx q[1];
rz(-0.60390632) q[1];
rz(-pi) q[2];
rz(-0.26594578) q[3];
sx q[3];
rz(-1.6308745) q[3];
sx q[3];
rz(1.657965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8043171) q[2];
sx q[2];
rz(-1.3108459) q[2];
sx q[2];
rz(2.6709225) q[2];
rz(-2.5110631) q[3];
sx q[3];
rz(-1.0258521) q[3];
sx q[3];
rz(1.7414909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077496342) q[0];
sx q[0];
rz(-0.12501669) q[0];
sx q[0];
rz(0.65573829) q[0];
rz(-0.31133044) q[1];
sx q[1];
rz(-0.89404023) q[1];
sx q[1];
rz(3.0701367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0229683) q[0];
sx q[0];
rz(-0.35425348) q[0];
sx q[0];
rz(-1.4262761) q[0];
rz(-1.2024443) q[2];
sx q[2];
rz(-1.2542274) q[2];
sx q[2];
rz(2.9376415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25117043) q[1];
sx q[1];
rz(-0.4931207) q[1];
sx q[1];
rz(1.432073) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3585594) q[3];
sx q[3];
rz(-2.0925412) q[3];
sx q[3];
rz(0.19475842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2283356) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(0.060001686) q[2];
rz(-1.2126728) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(-2.3771299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5948828) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(0.57719624) q[0];
rz(-2.3120841) q[1];
sx q[1];
rz(-1.0521768) q[1];
sx q[1];
rz(-0.83782354) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55623193) q[0];
sx q[0];
rz(-2.5149226) q[0];
sx q[0];
rz(0.13154948) q[0];
rz(-pi) q[1];
rz(-1.4805541) q[2];
sx q[2];
rz(-0.59207661) q[2];
sx q[2];
rz(-2.3088278) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8469593) q[1];
sx q[1];
rz(-2.3851352) q[1];
sx q[1];
rz(1.0119757) q[1];
rz(-pi) q[2];
rz(-2.2470705) q[3];
sx q[3];
rz(-1.2141879) q[3];
sx q[3];
rz(-1.5925533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.025658) q[2];
sx q[2];
rz(-1.5237153) q[2];
sx q[2];
rz(1.9177829) q[2];
rz(-0.4661679) q[3];
sx q[3];
rz(-0.40188447) q[3];
sx q[3];
rz(2.536072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.8849628) q[0];
sx q[0];
rz(-1.4361359) q[0];
sx q[0];
rz(-3.0404941) q[0];
rz(0.413232) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(1.0879999) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1032858) q[0];
sx q[0];
rz(-1.2438456) q[0];
sx q[0];
rz(-1.0590382) q[0];
rz(-pi) q[1];
rz(-1.9740943) q[2];
sx q[2];
rz(-0.77770644) q[2];
sx q[2];
rz(-2.4571927) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6551825) q[1];
sx q[1];
rz(-1.9014727) q[1];
sx q[1];
rz(2.8831456) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4721778) q[3];
sx q[3];
rz(-1.515404) q[3];
sx q[3];
rz(0.15633571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92855144) q[2];
sx q[2];
rz(-2.2553359) q[2];
sx q[2];
rz(1.586033) q[2];
rz(2.0689615) q[3];
sx q[3];
rz(-1.5242679) q[3];
sx q[3];
rz(0.13256375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9685386) q[0];
sx q[0];
rz(-0.88468495) q[0];
sx q[0];
rz(-1.7065077) q[0];
rz(2.3834719) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(3.039956) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1124737) q[0];
sx q[0];
rz(-2.2879507) q[0];
sx q[0];
rz(0.6458592) q[0];
rz(-pi) q[1];
rz(-2.0254214) q[2];
sx q[2];
rz(-1.2436507) q[2];
sx q[2];
rz(1.0516372) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5530049) q[1];
sx q[1];
rz(-1.8528474) q[1];
sx q[1];
rz(-1.4600919) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89438458) q[3];
sx q[3];
rz(-2.3830288) q[3];
sx q[3];
rz(2.0653084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.3397843) q[2];
sx q[2];
rz(-1.1016176) q[2];
sx q[2];
rz(-1.3327117) q[2];
rz(2.0594275) q[3];
sx q[3];
rz(-2.3887631) q[3];
sx q[3];
rz(1.4235206) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4448755) q[0];
sx q[0];
rz(-1.2514021) q[0];
sx q[0];
rz(-0.74309293) q[0];
rz(0.7630868) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(-0.9185763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615304) q[0];
sx q[0];
rz(-1.5788851) q[0];
sx q[0];
rz(1.5458406) q[0];
x q[1];
rz(0.95942504) q[2];
sx q[2];
rz(-2.583873) q[2];
sx q[2];
rz(3.0447247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6646361) q[1];
sx q[1];
rz(-1.0701706) q[1];
sx q[1];
rz(1.2229128) q[1];
rz(3.0094163) q[3];
sx q[3];
rz(-0.72849792) q[3];
sx q[3];
rz(2.8800509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4172198) q[2];
sx q[2];
rz(-1.1823187) q[2];
sx q[2];
rz(-2.7057538) q[2];
rz(-3.0271652) q[3];
sx q[3];
rz(-2.8947713) q[3];
sx q[3];
rz(2.54336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8082387) q[0];
sx q[0];
rz(-0.55753189) q[0];
sx q[0];
rz(0.58407855) q[0];
rz(2.7059879) q[1];
sx q[1];
rz(-1.6807115) q[1];
sx q[1];
rz(-1.5199419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3763894) q[0];
sx q[0];
rz(-2.0302773) q[0];
sx q[0];
rz(2.8110519) q[0];
rz(-pi) q[1];
rz(2.043739) q[2];
sx q[2];
rz(-0.45308896) q[2];
sx q[2];
rz(2.447213) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78306373) q[1];
sx q[1];
rz(-1.7768246) q[1];
sx q[1];
rz(1.5611783) q[1];
rz(-pi) q[2];
rz(-0.55840839) q[3];
sx q[3];
rz(-2.8388151) q[3];
sx q[3];
rz(-1.6807792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0674151) q[2];
sx q[2];
rz(-1.3520853) q[2];
sx q[2];
rz(1.1161067) q[2];
rz(-1.3793147) q[3];
sx q[3];
rz(-1.1032871) q[3];
sx q[3];
rz(0.11837676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24935687) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(-0.78999162) q[0];
rz(-0.063591592) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(-0.44073179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14549705) q[0];
sx q[0];
rz(-0.4410559) q[0];
sx q[0];
rz(-0.8121374) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67040261) q[2];
sx q[2];
rz(-0.60302654) q[2];
sx q[2];
rz(-0.5926026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4197582) q[1];
sx q[1];
rz(-1.4065308) q[1];
sx q[1];
rz(-2.052015) q[1];
x q[2];
rz(3.1001631) q[3];
sx q[3];
rz(-1.3815615) q[3];
sx q[3];
rz(1.7570868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8792087) q[2];
sx q[2];
rz(-2.3449506) q[2];
sx q[2];
rz(-2.4764496) q[2];
rz(-0.77196676) q[3];
sx q[3];
rz(-2.3699581) q[3];
sx q[3];
rz(1.5099248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36122286) q[0];
sx q[0];
rz(-0.64798111) q[0];
sx q[0];
rz(-1.004647) q[0];
rz(-1.7338344) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
sx q[1];
rz(-1.8515324) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0718075) q[0];
sx q[0];
rz(-1.8172853) q[0];
sx q[0];
rz(2.3355961) q[0];
x q[1];
rz(0.20310852) q[2];
sx q[2];
rz(-2.1611193) q[2];
sx q[2];
rz(3.1121569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5651144) q[1];
sx q[1];
rz(-1.2966411) q[1];
sx q[1];
rz(-1.3512011) q[1];
rz(-pi) q[2];
rz(-2.1716406) q[3];
sx q[3];
rz(-1.4267216) q[3];
sx q[3];
rz(0.010330095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1349692) q[2];
sx q[2];
rz(-1.9776191) q[2];
sx q[2];
rz(-0.25035614) q[2];
rz(2.1641459) q[3];
sx q[3];
rz(-1.2075295) q[3];
sx q[3];
rz(-1.6710963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866495) q[0];
sx q[0];
rz(-1.7087806) q[0];
sx q[0];
rz(1.9720672) q[0];
rz(0.74465887) q[1];
sx q[1];
rz(-0.79569334) q[1];
sx q[1];
rz(-2.4462499) q[1];
rz(-1.5834687) q[2];
sx q[2];
rz(-1.4274393) q[2];
sx q[2];
rz(1.9096331) q[2];
rz(-2.0674191) q[3];
sx q[3];
rz(-1.1840829) q[3];
sx q[3];
rz(0.89331762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
