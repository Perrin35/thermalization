OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(0.36261121) q[0];
rz(-2.2244722) q[1];
sx q[1];
rz(-2.6511104) q[1];
sx q[1];
rz(-2.7999556) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4877019) q[0];
sx q[0];
rz(-0.39429769) q[0];
sx q[0];
rz(2.8173692) q[0];
rz(-0.59092893) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(0.17609827) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4075549) q[1];
sx q[1];
rz(-1.2734379) q[1];
sx q[1];
rz(2.107891) q[1];
rz(-pi) q[2];
rz(1.1716236) q[3];
sx q[3];
rz(-1.7776383) q[3];
sx q[3];
rz(2.9575461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59387702) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(0.18125136) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618806) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(1.1454426) q[0];
rz(-1.068813) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(-0.72584814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94978588) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(1.9907065) q[0];
rz(2.0306573) q[2];
sx q[2];
rz(-2.1150555) q[2];
sx q[2];
rz(2.0366675) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1314288) q[1];
sx q[1];
rz(-1.5457075) q[1];
sx q[1];
rz(-1.2618834) q[1];
x q[2];
rz(-1.3594567) q[3];
sx q[3];
rz(-0.95892116) q[3];
sx q[3];
rz(-2.8733727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48250616) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-0.99728161) q[2];
rz(1.7287792) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(2.0130656) q[0];
rz(1.8380802) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(-2.2361141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826686) q[0];
sx q[0];
rz(-1.7610468) q[0];
sx q[0];
rz(-0.24567901) q[0];
rz(-1.1378023) q[2];
sx q[2];
rz(-2.275122) q[2];
sx q[2];
rz(-1.1906884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9943774) q[1];
sx q[1];
rz(-2.4265687) q[1];
sx q[1];
rz(-1.6836402) q[1];
rz(-pi) q[2];
rz(-3.0985673) q[3];
sx q[3];
rz(-2.283841) q[3];
sx q[3];
rz(-0.42792861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-2.1543489) q[2];
rz(0.045771249) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.1674081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2950738) q[0];
sx q[0];
rz(-0.82058883) q[0];
sx q[0];
rz(-1.2942765) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6378176) q[2];
sx q[2];
rz(-1.4132032) q[2];
sx q[2];
rz(2.4368844) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0427525) q[1];
sx q[1];
rz(-2.4272857) q[1];
sx q[1];
rz(-2.7318098) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79951841) q[3];
sx q[3];
rz(-0.82516731) q[3];
sx q[3];
rz(-2.094401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2771153) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(-2.5100822) q[2];
rz(-2.946092) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-2.6389129) q[0];
rz(-1.1133105) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(0.46708333) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.607476) q[0];
sx q[0];
rz(-0.42629888) q[0];
sx q[0];
rz(2.0712453) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.35891) q[2];
sx q[2];
rz(2.6168602) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.52006066) q[1];
sx q[1];
rz(-2.660775) q[1];
sx q[1];
rz(2.0562999) q[1];
rz(-pi) q[2];
rz(-0.15848666) q[3];
sx q[3];
rz(-2.2500258) q[3];
sx q[3];
rz(-2.2331626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-2.6494027) q[2];
rz(-0.28218937) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.870134) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(0.085993275) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(-1.496398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39280415) q[0];
sx q[0];
rz(-1.3747871) q[0];
sx q[0];
rz(2.1013573) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24106579) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(0.17503967) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3277668) q[1];
sx q[1];
rz(-1.5470978) q[1];
sx q[1];
rz(-1.2162672) q[1];
rz(0.71293998) q[3];
sx q[3];
rz(-1.5804277) q[3];
sx q[3];
rz(0.18129098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4795586) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(0.39624828) q[2];
rz(-0.59404343) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(-1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(2.0492045) q[0];
rz(-1.7842402) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(-0.011627442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41650018) q[0];
sx q[0];
rz(-1.9511127) q[0];
sx q[0];
rz(1.2854544) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3064012) q[2];
sx q[2];
rz(-3.1036658) q[2];
sx q[2];
rz(-1.7131625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36794084) q[1];
sx q[1];
rz(-0.546574) q[1];
sx q[1];
rz(-1.5252588) q[1];
rz(-pi) q[2];
x q[2];
rz(2.743268) q[3];
sx q[3];
rz(-2.4279804) q[3];
sx q[3];
rz(-2.6834965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45550436) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(-1.9912432) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31036723) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.9294552) q[0];
rz(-0.37777004) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6389248) q[0];
sx q[0];
rz(-1.272861) q[0];
sx q[0];
rz(1.7940679) q[0];
rz(0.059842589) q[2];
sx q[2];
rz(-1.4388196) q[2];
sx q[2];
rz(-3.0795385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2468176) q[1];
sx q[1];
rz(-2.2404039) q[1];
sx q[1];
rz(-1.0424022) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0231608) q[3];
sx q[3];
rz(-2.0698955) q[3];
sx q[3];
rz(-1.1008319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(1.2517694) q[2];
rz(2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-0.10929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733474) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(-2.9633203) q[0];
rz(-3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4561153) q[0];
sx q[0];
rz(-0.59887409) q[0];
sx q[0];
rz(-2.6045538) q[0];
rz(-pi) q[1];
rz(2.2237642) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(2.3248621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0563911) q[1];
sx q[1];
rz(-1.3136275) q[1];
sx q[1];
rz(2.609054) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2004553) q[3];
sx q[3];
rz(-0.59560532) q[3];
sx q[3];
rz(2.5411118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69958413) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(-2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853302) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(0.061696079) q[0];
rz(-1.8348947) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(1.0888938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058044) q[0];
sx q[0];
rz(-0.09011589) q[0];
sx q[0];
rz(-2.3401005) q[0];
rz(-0.18386545) q[2];
sx q[2];
rz(-0.92391787) q[2];
sx q[2];
rz(-1.136214) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1354462) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(1.0165434) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0763361) q[3];
sx q[3];
rz(-0.38443243) q[3];
sx q[3];
rz(-1.9264551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6178199) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(0.069996746) q[2];
rz(0.87674117) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(1.6246673) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(0.91296997) q[2];
sx q[2];
rz(-1.8987449) q[2];
sx q[2];
rz(1.434025) q[2];
rz(-1.2354479) q[3];
sx q[3];
rz(-1.609758) q[3];
sx q[3];
rz(2.5555425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
