OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(-2.3318113) q[0];
sx q[0];
rz(0.53139395) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(-1.2030503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42703585) q[0];
sx q[0];
rz(-1.7622951) q[0];
sx q[0];
rz(-1.3546076) q[0];
rz(2.4389624) q[2];
sx q[2];
rz(-2.812817) q[2];
sx q[2];
rz(-1.6871014) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0359218) q[1];
sx q[1];
rz(-2.8997313) q[1];
sx q[1];
rz(-0.24005228) q[1];
rz(-pi) q[2];
rz(0.12227998) q[3];
sx q[3];
rz(-2.8465726) q[3];
sx q[3];
rz(-2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8649341) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-2.9425088) q[2];
rz(1.52786) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(2.4075107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428225) q[0];
sx q[0];
rz(-0.12717371) q[0];
sx q[0];
rz(-0.81940991) q[0];
rz(-0.2858513) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(-1.8751289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3730777) q[0];
sx q[0];
rz(-0.74155945) q[0];
sx q[0];
rz(-0.57912796) q[0];
rz(2.0882294) q[2];
sx q[2];
rz(-0.45767637) q[2];
sx q[2];
rz(-2.6982754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4283735) q[1];
sx q[1];
rz(-1.1717147) q[1];
sx q[1];
rz(-2.748511) q[1];
x q[2];
rz(0.5760848) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9849898) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(-1.161969) q[2];
rz(0.087163838) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(-3.0676837) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53482985) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(2.8379922) q[0];
rz(-1.3820232) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(-2.4940925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8436463) q[0];
sx q[0];
rz(-0.93695153) q[0];
sx q[0];
rz(-1.514939) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0258255) q[2];
sx q[2];
rz(-2.2098944) q[2];
sx q[2];
rz(0.35312032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6806879) q[1];
sx q[1];
rz(-0.60084963) q[1];
sx q[1];
rz(-2.2520573) q[1];
rz(-0.59433523) q[3];
sx q[3];
rz(-0.96386516) q[3];
sx q[3];
rz(2.2567574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(1.5608609) q[2];
rz(-2.1598699) q[3];
sx q[3];
rz(-1.862062) q[3];
sx q[3];
rz(0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9075539) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(-2.2696944) q[0];
rz(2.7543228) q[1];
sx q[1];
rz(-2.498812) q[1];
sx q[1];
rz(0.29104582) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5440118) q[0];
sx q[0];
rz(-2.512092) q[0];
sx q[0];
rz(1.9660216) q[0];
rz(-pi) q[1];
rz(-2.8486738) q[2];
sx q[2];
rz(-1.2625164) q[2];
sx q[2];
rz(-2.3777547) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2070771) q[1];
sx q[1];
rz(-1.0672489) q[1];
sx q[1];
rz(0.7657004) q[1];
rz(2.1471359) q[3];
sx q[3];
rz(-1.7061702) q[3];
sx q[3];
rz(1.6614112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3691833) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(0.54405653) q[2];
rz(-2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-2.2756186) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13957025) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(-0.64506662) q[0];
rz(2.5158665) q[1];
sx q[1];
rz(-1.2236243) q[1];
sx q[1];
rz(-2.5114139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5298313) q[0];
sx q[0];
rz(-2.1453619) q[0];
sx q[0];
rz(-1.7163506) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1516045) q[2];
sx q[2];
rz(-0.79608166) q[2];
sx q[2];
rz(2.5615356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6781792) q[1];
sx q[1];
rz(-0.58600512) q[1];
sx q[1];
rz(1.9834118) q[1];
rz(-0.24124055) q[3];
sx q[3];
rz(-0.57096982) q[3];
sx q[3];
rz(2.8591448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3877635) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(-1.7774263) q[2];
rz(-1.2498614) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(-2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0548148) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-2.4247647) q[0];
rz(2.7038799) q[1];
sx q[1];
rz(-0.68044674) q[1];
sx q[1];
rz(0.81370083) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87156224) q[0];
sx q[0];
rz(-0.22261482) q[0];
sx q[0];
rz(1.2827669) q[0];
rz(-0.99102334) q[2];
sx q[2];
rz(-0.93587854) q[2];
sx q[2];
rz(0.53167508) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70255792) q[1];
sx q[1];
rz(-2.3137989) q[1];
sx q[1];
rz(2.031416) q[1];
rz(-pi) q[2];
rz(0.86088647) q[3];
sx q[3];
rz(-0.45300278) q[3];
sx q[3];
rz(-0.83356726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2563236) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.1876594) q[2];
rz(0.93196431) q[3];
sx q[3];
rz(-1.1340125) q[3];
sx q[3];
rz(-0.37117547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9880144) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(2.4555092) q[0];
rz(1.6185121) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(-2.4783321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30090573) q[0];
sx q[0];
rz(-1.5468883) q[0];
sx q[0];
rz(2.2211214) q[0];
rz(-1.3528321) q[2];
sx q[2];
rz(-2.1969165) q[2];
sx q[2];
rz(-1.8916212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38695947) q[1];
sx q[1];
rz(-0.52801758) q[1];
sx q[1];
rz(1.9792884) q[1];
x q[2];
rz(-2.5797964) q[3];
sx q[3];
rz(-1.1601163) q[3];
sx q[3];
rz(0.79515776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.474581) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(0.015080301) q[2];
rz(2.1298501) q[3];
sx q[3];
rz(-1.116131) q[3];
sx q[3];
rz(-0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-0.73260728) q[0];
sx q[0];
rz(-3.0342039) q[0];
rz(2.836851) q[1];
sx q[1];
rz(-1.3184897) q[1];
sx q[1];
rz(-0.4577786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37516884) q[0];
sx q[0];
rz(-2.1918115) q[0];
sx q[0];
rz(-2.7428438) q[0];
rz(-pi) q[1];
rz(-0.7465676) q[2];
sx q[2];
rz(-1.9779567) q[2];
sx q[2];
rz(0.50200576) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6785504) q[1];
sx q[1];
rz(-1.926534) q[1];
sx q[1];
rz(0.24865351) q[1];
rz(-pi) q[2];
rz(-2.5896766) q[3];
sx q[3];
rz(-1.2593927) q[3];
sx q[3];
rz(-0.50406721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7765939) q[2];
sx q[2];
rz(-1.4564617) q[2];
sx q[2];
rz(1.4314852) q[2];
rz(1.7163904) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9872221) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(-0.95348683) q[0];
rz(1.8158688) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(-0.58473933) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9033602) q[0];
sx q[0];
rz(-2.0582709) q[0];
sx q[0];
rz(1.4438629) q[0];
rz(-1.2051177) q[2];
sx q[2];
rz(-1.056991) q[2];
sx q[2];
rz(0.50525451) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.016420267) q[1];
sx q[1];
rz(-2.6965045) q[1];
sx q[1];
rz(2.9801286) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8381848) q[3];
sx q[3];
rz(-1.4328453) q[3];
sx q[3];
rz(1.7351868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3328302) q[2];
sx q[2];
rz(-0.82020438) q[2];
sx q[2];
rz(0.27627036) q[2];
rz(0.55109465) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(0.38366693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(-0.65548354) q[0];
rz(1.9845225) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(1.6190593) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009506) q[0];
sx q[0];
rz(-2.2146533) q[0];
sx q[0];
rz(-1.9615016) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4535975) q[2];
sx q[2];
rz(-1.3032459) q[2];
sx q[2];
rz(-2.5431395) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9138227) q[1];
sx q[1];
rz(-1.7350983) q[1];
sx q[1];
rz(2.5096202) q[1];
rz(-pi) q[2];
rz(-0.42201116) q[3];
sx q[3];
rz(-1.9415628) q[3];
sx q[3];
rz(-1.6251723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15554252) q[2];
sx q[2];
rz(-2.7337998) q[2];
sx q[2];
rz(-2.299451) q[2];
rz(-1.6048253) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(-3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.90606541) q[0];
sx q[0];
rz(-1.9530095) q[0];
sx q[0];
rz(-0.45146913) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(3.1149574) q[2];
sx q[2];
rz(-1.7799108) q[2];
sx q[2];
rz(-1.5540661) q[2];
rz(1.6668672) q[3];
sx q[3];
rz(-2.4484652) q[3];
sx q[3];
rz(-0.0948003) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
