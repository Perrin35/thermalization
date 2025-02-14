OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71275467) q[0];
sx q[0];
rz(-1.6102256) q[0];
sx q[0];
rz(1.1046326) q[0];
rz(-1.8072577) q[1];
sx q[1];
rz(5.5601064) q[1];
sx q[1];
rz(13.830166) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4975464) q[0];
sx q[0];
rz(-0.9263557) q[0];
sx q[0];
rz(0.30797663) q[0];
rz(-pi) q[1];
rz(-0.42277067) q[2];
sx q[2];
rz(-2.1029933) q[2];
sx q[2];
rz(-1.8025887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.059119) q[1];
sx q[1];
rz(-0.33519855) q[1];
sx q[1];
rz(-1.3223865) q[1];
rz(0.24407152) q[3];
sx q[3];
rz(-1.1435978) q[3];
sx q[3];
rz(2.0352767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9870712) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(1.7729574) q[2];
rz(2.9030419) q[3];
sx q[3];
rz(-2.7261901) q[3];
sx q[3];
rz(-3.0273048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40047613) q[0];
sx q[0];
rz(-1.1392765) q[0];
sx q[0];
rz(-1.9912632) q[0];
rz(-3.053983) q[1];
sx q[1];
rz(-1.8116415) q[1];
sx q[1];
rz(1.5706583) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99518665) q[0];
sx q[0];
rz(-1.0885462) q[0];
sx q[0];
rz(1.6907808) q[0];
x q[1];
rz(-2.803844) q[2];
sx q[2];
rz(-0.18924533) q[2];
sx q[2];
rz(0.35795975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24805574) q[1];
sx q[1];
rz(-0.94590294) q[1];
sx q[1];
rz(0.53433634) q[1];
rz(-pi) q[2];
rz(0.47111311) q[3];
sx q[3];
rz(-2.2918252) q[3];
sx q[3];
rz(-1.6278933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7972083) q[2];
sx q[2];
rz(-0.71343652) q[2];
sx q[2];
rz(0.1864645) q[2];
rz(-0.79948419) q[3];
sx q[3];
rz(-1.5461642) q[3];
sx q[3];
rz(2.1605087) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568273) q[0];
sx q[0];
rz(-1.3815877) q[0];
sx q[0];
rz(0.10936603) q[0];
rz(-0.60375396) q[1];
sx q[1];
rz(-1.8021288) q[1];
sx q[1];
rz(2.107479) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17149481) q[0];
sx q[0];
rz(-1.7037647) q[0];
sx q[0];
rz(1.7823969) q[0];
rz(-pi) q[1];
rz(1.6414406) q[2];
sx q[2];
rz(-1.6598668) q[2];
sx q[2];
rz(1.8235109) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2769985) q[1];
sx q[1];
rz(-1.2916451) q[1];
sx q[1];
rz(-1.6807081) q[1];
rz(-0.14520653) q[3];
sx q[3];
rz(-2.3646486) q[3];
sx q[3];
rz(2.8498385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5321396) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(0.54507059) q[2];
rz(-0.80101454) q[3];
sx q[3];
rz(-0.89371926) q[3];
sx q[3];
rz(3.0494704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48701778) q[0];
sx q[0];
rz(-1.6681404) q[0];
sx q[0];
rz(-0.45904485) q[0];
rz(1.1162988) q[1];
sx q[1];
rz(-0.79367677) q[1];
sx q[1];
rz(-0.8265411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3497884) q[0];
sx q[0];
rz(-1.1513745) q[0];
sx q[0];
rz(0.38695199) q[0];
x q[1];
rz(-3.0439348) q[2];
sx q[2];
rz(-1.6689166) q[2];
sx q[2];
rz(2.1922534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53163995) q[1];
sx q[1];
rz(-0.58216698) q[1];
sx q[1];
rz(-1.2083295) q[1];
x q[2];
rz(-1.6891805) q[3];
sx q[3];
rz(-1.2384129) q[3];
sx q[3];
rz(2.3207302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7894342) q[2];
sx q[2];
rz(-0.23098478) q[2];
sx q[2];
rz(-2.3918772) q[2];
rz(-1.1226783) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(-0.96021715) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66626755) q[0];
sx q[0];
rz(-0.98639494) q[0];
sx q[0];
rz(-3.0295897) q[0];
rz(0.22459596) q[1];
sx q[1];
rz(-1.1677531) q[1];
sx q[1];
rz(2.3602643) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3488779) q[0];
sx q[0];
rz(-1.166806) q[0];
sx q[0];
rz(-0.029026194) q[0];
rz(-pi) q[1];
rz(-0.26877578) q[2];
sx q[2];
rz(-0.34160638) q[2];
sx q[2];
rz(0.20526055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0722149) q[1];
sx q[1];
rz(-0.88000127) q[1];
sx q[1];
rz(1.0051954) q[1];
x q[2];
rz(2.1634401) q[3];
sx q[3];
rz(-1.7053805) q[3];
sx q[3];
rz(0.59813598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8319228) q[2];
sx q[2];
rz(-0.90193844) q[2];
sx q[2];
rz(2.208948) q[2];
rz(-2.0617088) q[3];
sx q[3];
rz(-0.44411689) q[3];
sx q[3];
rz(0.035695765) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90750736) q[0];
sx q[0];
rz(-1.1312753) q[0];
sx q[0];
rz(1.3908516) q[0];
rz(0.33175173) q[1];
sx q[1];
rz(-1.2650047) q[1];
sx q[1];
rz(-1.5914241) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54133139) q[0];
sx q[0];
rz(-2.2131767) q[0];
sx q[0];
rz(-3.0976901) q[0];
rz(-3.0671547) q[2];
sx q[2];
rz(-0.46069579) q[2];
sx q[2];
rz(1.6440533) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0968714) q[1];
sx q[1];
rz(-1.4145748) q[1];
sx q[1];
rz(-3.0169283) q[1];
rz(-pi) q[2];
rz(0.36339001) q[3];
sx q[3];
rz(-1.212032) q[3];
sx q[3];
rz(-0.10063688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3127689) q[2];
sx q[2];
rz(-0.88461107) q[2];
sx q[2];
rz(-1.5213607) q[2];
rz(-0.96794266) q[3];
sx q[3];
rz(-1.5588375) q[3];
sx q[3];
rz(2.2255285) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5172326) q[0];
sx q[0];
rz(-2.7150798) q[0];
sx q[0];
rz(2.970001) q[0];
rz(1.0393556) q[1];
sx q[1];
rz(-1.2587222) q[1];
sx q[1];
rz(-1.7003869) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489451) q[0];
sx q[0];
rz(-1.2496557) q[0];
sx q[0];
rz(-1.7016861) q[0];
rz(-0.19564678) q[2];
sx q[2];
rz(-0.52426978) q[2];
sx q[2];
rz(-1.9877246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4644147) q[1];
sx q[1];
rz(-0.59301335) q[1];
sx q[1];
rz(-0.89927425) q[1];
x q[2];
rz(-1.3492382) q[3];
sx q[3];
rz(-2.5668813) q[3];
sx q[3];
rz(1.2971085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8695996) q[2];
sx q[2];
rz(-1.7337493) q[2];
sx q[2];
rz(0.54455152) q[2];
rz(-2.6416687) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(2.669529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8498103) q[0];
sx q[0];
rz(-2.6047459) q[0];
sx q[0];
rz(0.52919069) q[0];
rz(-2.5657907) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(-1.1189438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3258695) q[0];
sx q[0];
rz(-1.6224321) q[0];
sx q[0];
rz(1.5918094) q[0];
x q[1];
rz(0.37092692) q[2];
sx q[2];
rz(-1.7417522) q[2];
sx q[2];
rz(1.8425187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0487489) q[1];
sx q[1];
rz(-2.1511937) q[1];
sx q[1];
rz(1.012085) q[1];
x q[2];
rz(2.9319134) q[3];
sx q[3];
rz(-2.0884313) q[3];
sx q[3];
rz(0.28214327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.240856) q[2];
sx q[2];
rz(-2.6764328) q[2];
sx q[2];
rz(-0.71211234) q[2];
rz(1.6093048) q[3];
sx q[3];
rz(-2.1963547) q[3];
sx q[3];
rz(-1.7942662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.79494548) q[0];
sx q[0];
rz(-2.0136588) q[0];
sx q[0];
rz(1.0711063) q[0];
rz(1.2062997) q[1];
sx q[1];
rz(-2.1251528) q[1];
sx q[1];
rz(1.4869022) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95222118) q[0];
sx q[0];
rz(-2.4444488) q[0];
sx q[0];
rz(2.8482262) q[0];
rz(-pi) q[1];
rz(0.63196147) q[2];
sx q[2];
rz(-2.2990531) q[2];
sx q[2];
rz(1.2974844) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35569977) q[1];
sx q[1];
rz(-1.744387) q[1];
sx q[1];
rz(-0.29533556) q[1];
rz(-2.0196718) q[3];
sx q[3];
rz(-1.316183) q[3];
sx q[3];
rz(-1.5625579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6335166) q[2];
sx q[2];
rz(-0.8997007) q[2];
sx q[2];
rz(3.0637975) q[2];
rz(0.22732321) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(-2.0859065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530497) q[0];
sx q[0];
rz(-0.74497861) q[0];
sx q[0];
rz(2.7689834) q[0];
rz(0.44652069) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(2.2850697) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4093709) q[0];
sx q[0];
rz(-1.5268832) q[0];
sx q[0];
rz(-3.0847021) q[0];
rz(-3.0356016) q[2];
sx q[2];
rz(-1.4376831) q[2];
sx q[2];
rz(3.1337381) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7584655) q[1];
sx q[1];
rz(-1.9680319) q[1];
sx q[1];
rz(-3.0137311) q[1];
rz(2.4803745) q[3];
sx q[3];
rz(-2.2228129) q[3];
sx q[3];
rz(-2.2673502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7591758) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(0.55553931) q[2];
rz(-2.7555452) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(0.66421318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0059218) q[0];
sx q[0];
rz(-1.6914524) q[0];
sx q[0];
rz(1.2941262) q[0];
rz(0.80264965) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(-0.35292179) q[2];
sx q[2];
rz(-1.6833932) q[2];
sx q[2];
rz(2.0298454) q[2];
rz(3.0714645) q[3];
sx q[3];
rz(-0.84134103) q[3];
sx q[3];
rz(0.76920912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
