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
rz(2.6142081) q[0];
sx q[0];
rz(-2.5684147) q[0];
sx q[0];
rz(2.5249124) q[0];
rz(-1.0083899) q[1];
sx q[1];
rz(-2.402306) q[1];
sx q[1];
rz(-0.33831212) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36114281) q[0];
sx q[0];
rz(-2.8964213) q[0];
sx q[0];
rz(0.14921363) q[0];
rz(2.7213547) q[2];
sx q[2];
rz(-0.8367238) q[2];
sx q[2];
rz(-0.24009304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.19060082) q[1];
sx q[1];
rz(-1.436232) q[1];
sx q[1];
rz(-2.5959204) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81387384) q[3];
sx q[3];
rz(-1.699109) q[3];
sx q[3];
rz(1.8112884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3931291) q[2];
sx q[2];
rz(-1.4145565) q[2];
sx q[2];
rz(0.65790042) q[2];
rz(-0.45514485) q[3];
sx q[3];
rz(-0.15076605) q[3];
sx q[3];
rz(-1.8337102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.2748134) q[0];
sx q[0];
rz(-1.7300737) q[0];
sx q[0];
rz(0.28208062) q[0];
rz(2.8981949) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(2.3928221) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215749) q[0];
sx q[0];
rz(-0.87792464) q[0];
sx q[0];
rz(1.3189032) q[0];
rz(-pi) q[1];
rz(3.1339945) q[2];
sx q[2];
rz(-1.4491334) q[2];
sx q[2];
rz(-1.0341782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2255324) q[1];
sx q[1];
rz(-2.6301363) q[1];
sx q[1];
rz(-1.3813853) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5010819) q[3];
sx q[3];
rz(-1.3718714) q[3];
sx q[3];
rz(0.45123842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8613209) q[2];
sx q[2];
rz(-1.6161641) q[2];
sx q[2];
rz(-1.8005499) q[2];
rz(-1.1567814) q[3];
sx q[3];
rz(-0.78514922) q[3];
sx q[3];
rz(-2.3780499) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083104) q[0];
sx q[0];
rz(-2.9756727) q[0];
sx q[0];
rz(2.4842343) q[0];
rz(-1.9961458) q[1];
sx q[1];
rz(-0.95445389) q[1];
sx q[1];
rz(-0.81833902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2812202) q[0];
sx q[0];
rz(-1.91433) q[0];
sx q[0];
rz(1.8247129) q[0];
rz(1.7162227) q[2];
sx q[2];
rz(-2.8231695) q[2];
sx q[2];
rz(2.2216036) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9200478) q[1];
sx q[1];
rz(-2.2755879) q[1];
sx q[1];
rz(-0.66001604) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0243689) q[3];
sx q[3];
rz(-1.6370341) q[3];
sx q[3];
rz(0.57285786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98693097) q[2];
sx q[2];
rz(-1.8795452) q[2];
sx q[2];
rz(1.7884802) q[2];
rz(-0.96495572) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(2.7195948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44593909) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(1.9406142) q[0];
rz(-3.1383842) q[1];
sx q[1];
rz(-1.708958) q[1];
sx q[1];
rz(-0.23689717) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0081404765) q[0];
sx q[0];
rz(-2.4137375) q[0];
sx q[0];
rz(-0.86489622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9658674) q[2];
sx q[2];
rz(-0.64470664) q[2];
sx q[2];
rz(2.3852661) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8767406) q[1];
sx q[1];
rz(-0.63794604) q[1];
sx q[1];
rz(-1.2400342) q[1];
rz(-pi) q[2];
rz(-0.75613344) q[3];
sx q[3];
rz(-2.2069158) q[3];
sx q[3];
rz(2.3622733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0420456) q[2];
sx q[2];
rz(-1.2914265) q[2];
sx q[2];
rz(2.715204) q[2];
rz(2.1060627) q[3];
sx q[3];
rz(-1.6798881) q[3];
sx q[3];
rz(-2.5199913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7972888) q[0];
sx q[0];
rz(-0.042003691) q[0];
sx q[0];
rz(-2.7916743) q[0];
rz(-1.2087076) q[1];
sx q[1];
rz(-1.8588926) q[1];
sx q[1];
rz(0.0016317687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8953319) q[0];
sx q[0];
rz(-0.29731942) q[0];
sx q[0];
rz(0.60054442) q[0];
rz(1.5320184) q[2];
sx q[2];
rz(-1.2735521) q[2];
sx q[2];
rz(0.69678885) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.220158) q[1];
sx q[1];
rz(-1.5234656) q[1];
sx q[1];
rz(-1.7478075) q[1];
rz(-pi) q[2];
rz(2.9309421) q[3];
sx q[3];
rz(-1.306476) q[3];
sx q[3];
rz(1.5624865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1899167) q[2];
sx q[2];
rz(-0.92477551) q[2];
sx q[2];
rz(-2.9803989) q[2];
rz(2.7023884) q[3];
sx q[3];
rz(-2.3930211) q[3];
sx q[3];
rz(0.85328931) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8301903) q[0];
sx q[0];
rz(-1.4051733) q[0];
sx q[0];
rz(2.3468974) q[0];
rz(1.7757802) q[1];
sx q[1];
rz(-2.3410485) q[1];
sx q[1];
rz(-0.47439233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649523) q[0];
sx q[0];
rz(-0.65529233) q[0];
sx q[0];
rz(2.220015) q[0];
x q[1];
rz(-1.2016868) q[2];
sx q[2];
rz(-1.0449191) q[2];
sx q[2];
rz(-0.78247082) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77955708) q[1];
sx q[1];
rz(-2.1367225) q[1];
sx q[1];
rz(1.0841682) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0211518) q[3];
sx q[3];
rz(-2.4826038) q[3];
sx q[3];
rz(3.0699108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5279493) q[2];
sx q[2];
rz(-1.2443292) q[2];
sx q[2];
rz(2.6141686) q[2];
rz(-2.534965) q[3];
sx q[3];
rz(-2.9546723) q[3];
sx q[3];
rz(2.0691779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793176) q[0];
sx q[0];
rz(-2.9476705) q[0];
sx q[0];
rz(0.79757565) q[0];
rz(-0.89353117) q[1];
sx q[1];
rz(-0.83502665) q[1];
sx q[1];
rz(-0.28018793) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55927568) q[0];
sx q[0];
rz(-1.7452551) q[0];
sx q[0];
rz(1.972439) q[0];
rz(2.3827219) q[2];
sx q[2];
rz(-0.810383) q[2];
sx q[2];
rz(-2.9019865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57443212) q[1];
sx q[1];
rz(-0.48266294) q[1];
sx q[1];
rz(-1.5437276) q[1];
rz(-pi) q[2];
rz(0.97902043) q[3];
sx q[3];
rz(-1.8734402) q[3];
sx q[3];
rz(-1.7906063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7402652) q[2];
sx q[2];
rz(-1.836931) q[2];
sx q[2];
rz(-1.2124445) q[2];
rz(-2.9017743) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(0.56813204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3577393) q[0];
sx q[0];
rz(-1.4171866) q[0];
sx q[0];
rz(1.3215815) q[0];
rz(0.18121885) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(-0.063057335) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94751271) q[0];
sx q[0];
rz(-2.3499858) q[0];
sx q[0];
rz(0.20157728) q[0];
rz(-pi) q[1];
rz(-1.8562154) q[2];
sx q[2];
rz(-1.6550502) q[2];
sx q[2];
rz(0.38049305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93062295) q[1];
sx q[1];
rz(-1.9038868) q[1];
sx q[1];
rz(-1.3107783) q[1];
x q[2];
rz(0.0029011528) q[3];
sx q[3];
rz(-1.408268) q[3];
sx q[3];
rz(1.8123466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41398373) q[2];
sx q[2];
rz(-1.0980282) q[2];
sx q[2];
rz(-1.46924) q[2];
rz(-1.3937048) q[3];
sx q[3];
rz(-1.7017476) q[3];
sx q[3];
rz(2.9318455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53836981) q[0];
sx q[0];
rz(-2.5312238) q[0];
sx q[0];
rz(2.1446153) q[0];
rz(-0.793055) q[1];
sx q[1];
rz(-1.4184364) q[1];
sx q[1];
rz(2.0319895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0080896688) q[0];
sx q[0];
rz(-2.4937196) q[0];
sx q[0];
rz(-2.6459207) q[0];
rz(-pi) q[1];
rz(-1.1446321) q[2];
sx q[2];
rz(-2.2467016) q[2];
sx q[2];
rz(0.083783178) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0834137) q[1];
sx q[1];
rz(-1.021487) q[1];
sx q[1];
rz(0.55725248) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.505198) q[3];
sx q[3];
rz(-1.706746) q[3];
sx q[3];
rz(3.1042447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41254607) q[2];
sx q[2];
rz(-1.46526) q[2];
sx q[2];
rz(-1.868978) q[2];
rz(-2.0160969) q[3];
sx q[3];
rz(-0.19386217) q[3];
sx q[3];
rz(3.061749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37128714) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(3.0433997) q[0];
rz(-1.498361) q[1];
sx q[1];
rz(-1.9941092) q[1];
sx q[1];
rz(-2.3032545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1824552) q[0];
sx q[0];
rz(-2.3172624) q[0];
sx q[0];
rz(-1.3289544) q[0];
x q[1];
rz(-1.5151843) q[2];
sx q[2];
rz(-1.8225553) q[2];
sx q[2];
rz(2.1625569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0948719) q[1];
sx q[1];
rz(-1.1390242) q[1];
sx q[1];
rz(-0.31109613) q[1];
rz(-pi) q[2];
rz(2.7097852) q[3];
sx q[3];
rz(-1.7718959) q[3];
sx q[3];
rz(2.3922298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3457634) q[2];
sx q[2];
rz(-0.886262) q[2];
sx q[2];
rz(0.32540992) q[2];
rz(0.77643967) q[3];
sx q[3];
rz(-1.0151981) q[3];
sx q[3];
rz(2.5344892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14071874) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(2.8554032) q[1];
sx q[1];
rz(-1.9356526) q[1];
sx q[1];
rz(-2.0331358) q[1];
rz(-0.76079615) q[2];
sx q[2];
rz(-1.2789055) q[2];
sx q[2];
rz(0.61550752) q[2];
rz(-0.48578942) q[3];
sx q[3];
rz(-2.494907) q[3];
sx q[3];
rz(-3.1068556) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
